#!/usr/bin/env python3
import json
import os
import sys
from typing import Literal

import requests
from google import genai
from google.genai import types
from pydantic import BaseModel, Field

# Configurations from environment
GITHUB_TOKEN = os.getenv("GITHUB_TOKEN")
GEMINI_API_KEY = os.getenv("GEMINI_API_KEY")
REPO = os.getenv("GITHUB_REPOSITORY")
PR_NUMBER = os.getenv("PR_NUMBER")
MODEL_NAME = "gemini-3.8-flash"

SEVERITY_BADGE = {
    "blocking": "🔴 **BLOCKING**",
    "major": "🟠 **MAJOR**",
    "minor": "🟡 **MINOR**",
}

if not all([GITHUB_TOKEN, GEMINI_API_KEY, REPO, PR_NUMBER]):
    raise ValueError("Missing required environment variables.")


def get_pr_metadata():
    """Fetches the PR title and description body."""
    url = f"https://api.github.com/repos/{REPO}/pulls/{PR_NUMBER}"
    headers = {
        "Authorization": f"token {GITHUB_TOKEN}",
        "Accept": "application/vnd.github.v3+json",
    }
    response = requests.get(url, headers=headers)
    response.raise_for_status()
    data = response.json()
    return {"title": data.get("title", ""), "body": data.get("body", "")}


def get_pr_diff():
    response = requests.get(
        f"https://patch-diff.githubusercontent.com/raw/{REPO}/pull/{PR_NUMBER}.diff"
    )
    response.raise_for_status()
    return response.text


def post_review_comments(risk_summary, comments):
    url = f"https://api.github.com/repos/{REPO}/pulls/{PR_NUMBER}/reviews"
    headers = {
        "Authorization": f"token {GITHUB_TOKEN}",
        "Accept": "application/vnd.github.v3+json",
    }

    github_comments = []
    has_blocking = False

    # Sort comments in case one fails.
    severity_order = {"blocking": 0, "major": 1, "minor": 2}
    comments.sort(key=lambda c: severity_order.get(c.get("severity", "minor"), 3))

    for c in comments:
        severity = c.get("severity", "minor")
        has_blocking = severity == "blocking"

        badge = SEVERITY_BADGE.get(severity, "")
        body = f"{badge}\n\n{c['explanation']}" if badge else c["explanation"]

        suggestion = c.get("suggestion")
        if suggestion:
            body += f"\n\n```suggestion\n{suggestion.strip()}\n```"

        github_comments.append(
            {
                "path": c["path"],
                "line": int(c["line"]),
                "side": "RIGHT",
                "body": body,
            }
        )

    # Escalate the review event if anything blocking was found
    event = "REQUEST_CHANGES" if has_blocking else "COMMENT"

    payload = {
        "body": f"# 🤖 Automated Gemini Code Review\n\n{risk_summary}",
        "event": event,
        "comments": github_comments,
    }

    res = requests.post(url, headers=headers, json=payload)
    res.raise_for_status()
    print(f"Successfully posted {len(github_comments)} review comments ({event}).")


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
    severity: Literal["blocking", "major", "minor"] = Field(
        description="blocking = would cause incorrect physics/results or a crash; "
        "major = correctness risk needing verification; minor = efficiency/robustness"
    )
    suggestion: str | None = Field(
        default=None,
        description="Exact code replacement, ONLY if you are confident in a concrete fix. "
        "Omit rather than fabricate a plausible-looking but unverified fix.",
    )


class ReviewPayload(BaseModel):
    risk_summary: str = Field(
        description="2-4 sentences: what subsystems this PR touches and the overall risk level. Do not restate the PR description — assume the reader has already read it."
    )
    comments: list[ReviewComment] = Field(description="List of inline code suggestions")


SYSTEM_INSTRUCTION = """
You are providing a strict and uncompromising pull request code review for Nyx, a high-fidelity, fast, and validated astrodynamics toolkit.
The toolkit is written in Rust with Python bindings via PyO3.

Your knowledge of astrodynamics, mission design and orbit determination spans from encyclopedias to
the latest state-of-the-art findings from AIAA/AAS Astrodynamics Specialist Conference papers.
You can accurately reference specific sections of the JPL DESCANSO monographs, the Ansys STK documentation, the ODTK MathSpec,
the CCSDS Blue Books, and other references of similar caliber. When flagging a physics/algorithm concern, name the specific reference and section that justifies the correct approach (e.g., 'Vallado 4th ed. §3.7' or 'DESCANSO Monograph 8, Ch. 4').

Nyx uses ANISE for all SPICE-related computation, frame transformations, rotation calculations, orbital element calculations, etc.
ANISE is a thread-safe, zero-cost alternative to NASA SPICE toolkit, computing spacecraft, planetary, coordinate frame, instrument transformations,
ground station visibility, and orbital elements. Engineered for high-throughput Python concurrent execution and flight software,
with proven lunar flight heritage on the Firefly Blue Ghost lunar lander, and on a large scale Earth orbiting constellation.
All time scale and duration computations are managed through Hifitime. Hifitime is an overflow-safe, high-performance datetime library providing
leap-second-correct nanosecond precision across UTC, GPST, and relativistic time-scales. Flight-proven in lunar and deep-space missions.

CRITICAL METRIC: You are evaluated solely on identifying architectural flaws and physical/mathematical errors or shortcuts that may have catastrophic effects in flight.
STRICTLY FORBIDDEN: Do not comment on code formatting, style variations, documentation formatting, or trivial typos in comments. If a change does not risk breaking execution, thread safety, serialization, or physical precision, IGNORE IT.

It is a correct and expected outcome to return zero comments for a file with no physics, serialization, or interface issues. Precision matters more than volume — a false positive that blocks a valid PR is worse than a missed nitpick. Do not manufacture a finding to appear thorough.

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

    # The Interactions API flattens configuration parameters and unifies model communication
    interaction = client.interactions.create(
        model=MODEL_NAME,
        input=prompt_content,
        system_instruction=SYSTEM_INSTRUCTION,
        response_format={
            "type": "text",
            "mime_type": "application/json",
            "schema": ReviewPayload,
        },
        generation_config={"temperature": 0.5, "thinking_level": "high"},
    )

    try:
        # Access the text output directly from interaction.output_text
        review_data = json.loads(interaction.output_text)
        post_review_comments(
            review_data.get("risk_summary", "No summary provided"),
            review_data.get("comments", []),
        )
    except (json.JSONDecodeError, AttributeError):
        print("Failed to parse model output or post review.")
        print(
            f"Raw response: {getattr(interaction, 'output_text', 'No output text found')}"
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
