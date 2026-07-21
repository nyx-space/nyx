#!/usr/bin/env python3
import os
import sys
import json
import requests
from google import genai
from google.genai import types
from pydantic import BaseModel, Field

# Configurations from environment
GITHUB_TOKEN = os.getenv("GITHUB_TOKEN")
GEMINI_API_KEY = os.getenv("GEMINI_API_KEY")
REPO = os.getenv("GITHUB_REPOSITORY")
PR_NUMBER = os.getenv("PR_NUMBER")
MODEL_NAME = "gemini-3.1-pro-preview"

if not all([GITHUB_TOKEN, GEMINI_API_KEY, REPO, PR_NUMBER]):
    raise ValueError("Missing required environment variables.")


def get_pr_metadata():
    """Fetches the PR title and description body."""
    url = f"https://api.github.com/repos/{REPO}/pulls/{PR_NUMBER}"
    headers = {
        "Authorization": f"token {GITHUB_TOKEN}",
        "Accept": "application.vnd.github.v3+json",
    }
    response = requests.get(url, headers=headers)
    response.raise_for_status()
    data = response.json()
    return {"title": data.get("title", ""), "body": data.get("body", "")}


def get_pr_diff():
    response = requests.get(f"https://patch-diff.githubusercontent.com/raw/nyx-space/nyx/pull/{PR_NUMBER}.patch")
    response.raise_for_status()
    return response.text


def post_review_comments(summary_markdown, comments):
    url = f"https://api.github.com/repos/{REPO}/pulls/{PR_NUMBER}/reviews"
    headers = {
        "Authorization": f"token {GITHUB_TOKEN}",
        "Accept": "application.vnd.github.v3+json",
    }

    # Map the JSON structure to GitHub's Review Comments API format
    github_comments = []
    for c in comments:
        github_comments.append(
            {
                "path": c["path"],
                "line": int(c["line"]),
                "side": "RIGHT",
                "body": f"{c['explanation']}\n\n```suggestion\n{c['suggestion'].strip()}\n```",
            }
        )

    payload = {
        "body": f"# 🤖 Automated Gemini Code Review\n\n{summary_markdown}",
        "event": "COMMENT",
        "comments": github_comments,
    }

    res = requests.post(url, headers=headers, json=payload)
    res.raise_for_status()
    print(f"Successfully posted {len(github_comments)} review comments.")


# Define standard Pydantic models for structured output
class ReviewComment(BaseModel):
    path: str = Field(
        description="Relative file path from repo root. Must match the exact file header in the diff."
    )
    line: int = Field(
        description="The exact target line number in the NEW version of the file where the modification occurs. Crucial: Must be a line present within the provided diff hunk lines."
    )
    explanation: str = Field(
        description="Concise architectural rationale for the change."
    )
    suggestion: str = Field(
        description="The exact code replacement block. Do not include markdown wrappers here."
    )

class ReviewPayload(BaseModel):
    summary_markdown: str = Field(
        description="The high-level PR summary formatted exactly according to the requested markdown template."
    )
    comments: list[ReviewComment] = Field(description="List of inline code suggestions")

SYSTEM_INSTRUCTION = """
You are providing a strict and uncompromising pull request code review for Nyx, a high-fidelity, fast, and validated astrodynamics toolkit.
The toolkit is written in Rust with Python bindings via PyO3.

Your knowledge of astrodynamics, mission design and orbit determination spans from encyclopedias to
the latest state-of-the-art findings from AIAA/AAS Astrodynamics Specialist Conference papers.
You can accurately reference specific sections of the JPL DESCANSO monographs, the Ansys STK documentation, the ODTK MathSpec,
the CCSDS Blue Books, and other references of similar caliber.

Nyx uses ANISE for all SPICE-related computation, frame transformations, rotation calculations, orbital element calculations, etc.
ANISE is a thread-safe, zero-cost alternative to NASA SPICE toolkit, computing spacecraft, planetary, coordinate frame, instrument transformations,
ground station visibility, and orbital elements. Engineered for high-throughput Python concurrent execution and flight software,
with proven lunar flight heritage on the Firefly Blue Ghost lunar lander, and on a large scale Earth orbiting constellation.
All time scale and duration computations are managed through Hifitime. Hifitime is an overflow-safe, high-performance datetime library providing
leap-second-correct nanosecond precision across UTC, GPST, and relativistic time-scales. Flight-proven in lunar and deep-space missions.

CRITICAL METRIC: You are evaluated solely on identifying architectural flaws and physical/mathematical errors or shortcuts that may have catastrophic effects in flight.
STRICTLY FORBIDDEN: Do not comment on code formatting, style variations, documentation formatting, or trivial typos in comments. If a change does not risk breaking execution, thread safety, serialization, or physical precision, IGNORE IT.

Target Evaluation Priorities:
1. Astrodynamics & Physics: Evaluate the underlying math and physics. Cross-reference implementation against state-of-the-art methods. Flag invalid assumptions, dangerous simplifications, or numerical instability risks.
2. Serialization Integrity: When a data structure is altered, trace its initialization and serialization footprints. Look for manual implementations of configuration traits. Flag missing field map additions causing schema mismatches.
3. Cross-Language Interface Invariants (Rust/Python): Scrutinize boundaries where native Rust logic meets Python bindings (`#[pyclass]`, `#[pymethods]`). Verify type compatibility.
4. Memory Allocations & Efficiency: Flag non-zero-cost abstractions, explicit vector allocations, boxes, or `.clone()` invocations inside tight numerical propagation loops.
5. Chronometry & Kinematics: Enforce rigorous verification of time scale conversions using hifitime and coordinate systems via ANISE.

You must output valid JSON matching the schema precisely. Do not hallucinate line numbers. If a file requires no changes, omit it from the array.
You will be evaluated based on the absolute structural accuracy of your line targets.

CRITICAL LINE-NUMBER DIRECTIVE:
- For every review comment you generate, the `line` property MUST correspond strictly to a valid line number added or modified in the NEW file context as presented in the unified diff headers (`@@ -... +... @@`).
- If a structural omission occurs (e.g., a field was omitted from an array or trait map downstream in the file), place the recommendation directly on the closest modification line or instantiation block visible within that specific diff hunk. Never target lines outside the provided hunks.

TEMPLATE INSTRUCTION:
You must also generate a high-level summary of the pull request using EXACTLY the following Markdown template. Fill in the sections appropriately based on the diff. If a section has no changes, leave it as "No change" (or the default text provided in the template). Do not alter the headings.

```markdown
## Summary

**Summarize the proposed changes**

## Architectural Changes

<!-- List any architectural changes made in this pull request, including any changes to the directory structure, file organization, or dependencies. -->

No change

## New Features

<!-- List any new features added in this pull request, including any new tools or functionality. -->

No change

## Improvements

<!-- List any improvements made in this pull request, including any performance optimizations, bug fixes, or other enhancements. -->

No change

## Bug Fixes

<!-- List any bug fixes made in this pull request, including any issues that were resolved. -->

No change

## Testing and validation

<!-- Please provide information on how the changes in this pull request were tested, including any new tests that were added or existing tests that were modified. -->

**Detail the changes in tests, including new tests and validations**

## Documentation

<!-- Detail documentation changes if this pull request primarily deals with documentation. -->

This PR does not primarily deal with documentation changes.
```
"""

def main():
    metadata = get_pr_metadata()
    diff_data = get_pr_diff()

    # Bundle description intent with the code changes
    prompt_content = f"""
=== PULL REQUEST METADATA ===
Title: {metadata["title"]}
Description:
{metadata["body"]}

=== UNIFIED CODE DIFF ===
{diff_data}
"""

    client = genai.Client(api_key=GEMINI_API_KEY)

    response = client.models.generate_content(
        model=MODEL_NAME,
        contents=prompt_content,
        config=types.GenerateContentConfig(
            system_instruction=SYSTEM_INSTRUCTION,
            response_mime_type="application/json",
            response_schema=ReviewPayload,
            temperature=0.1,  # Low temperature minimizes creative hallucination of lines
        ),
    )

    try:
        review_data = json.loads(response.text)
        post_review_comments(review_data.get("comments", []))
    except Exception as e:
        print(f"Failed to parse model output or post review: {e}")
        print(f"Raw response: {response.text}")
        sys.exit(1)


if __name__ == "__main__":
    main()
