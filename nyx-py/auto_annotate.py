import os
import re
import subprocess


def rust_to_py_type(rust_type: str, class_name: str) -> str:
    rust_type = rust_type.strip()
    if not rust_type:
        return "None"

    # Handle Result types: Result<T, E> or PyResult<T>
    if rust_type.startswith("Result<") or rust_type.startswith("PyResult<"):
        depth = 0
        inner = ""
        start_idx = rust_type.index("<") + 1
        for i in range(start_idx, len(rust_type)):
            c = rust_type[i]
            if c == "<":
                depth += 1
            elif c == ">":
                if depth == 0:
                    break
                depth -= 1
            elif c == "," and depth == 0:
                break
            inner += c
        rust_type = inner.strip()

    # Handle Option types
    if rust_type.startswith("Option<"):
        inner = rust_type[7:-1].strip()
        py_inner = rust_to_py_type(inner, class_name)
        if py_inner == "None":
            return "None"
        return f"{py_inner} | None"

    # Handle Vec / arrays
    if rust_type.startswith("Vec<") or rust_type.startswith("&["):
        if rust_type.startswith("Vec<"):
            inner = rust_type[4:-1].strip()
        else:
            inner = rust_type[2:-1].strip()
        py_inner = rust_to_py_type(inner, class_name)
        return f"list[{py_inner}]"

    # Stripping references, Lifetimes, mut, etc.
    rust_type = (
        rust_type.replace("&mut ", "")
        .replace("&self", "")
        .replace("&mut self", "")
        .replace("&", "")
        .strip()
    )
    if "'" in rust_type:
        rust_type = (
            re.sub(r"'[a-zA-Z_0-9]+", "", rust_type)
            .replace("< ,", "<")
            .replace("<,", "<")
            .strip()
        )

    if rust_type == "Self" or rust_type == "self":
        if class_name.startswith("Py") and class_name != "PyClass":
            return class_name[2:]
        return class_name

    if rust_type.startswith("Bound<"):
        m = re.match(r"Bound\s*<\s*(?:[^,]+,\s*)?([a-zA-Z0-9_<>\s]+)\s*>", rust_type)
        if m:
            rust_type = m.group(1).strip()
        else:
            return "typing.Any"

    # Standard mappings
    mappings = {
        "String": "str",
        "str": "str",
        "f64": "float",
        "f32": "float",
        "usize": "int",
        "u64": "int",
        "i64": "int",
        "u32": "int",
        "i32": "int",
        "u16": "int",
        "i16": "int",
        "u8": "int",
        "i8": "int",
        "bool": "bool",
        "()": "None",
        "None": "None",
    }

    if rust_type in mappings:
        return mappings[rust_type]

    if rust_type.startswith("Py") and len(rust_type) > 2 and rust_type[2].isupper():
        return rust_type[2:]

    return rust_type


def get_missing_rtypes(module_name: str) -> list[tuple[str, str]]:
    cmd = [".venv/bin/python", "generate_stubs.py", module_name, "temp.pyi"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    warnings = []
    # Search in both stdout and stderr
    combined_output = res.stdout + "\n" + res.stderr
    for line in combined_output.splitlines():
        # format: "Warning: The return type of nyx_space.mission_design.AtmDensity.earth_exponential is missing from the function documentation"
        if "The return type of" in line and "is missing" in line:
            parts = (
                line.split("The return type of ")[1].split(" is missing")[0].split(".")
            )
            if len(parts) >= 4:
                class_name = parts[-2]
                method_name = parts[-1]
                warnings.append((class_name, method_name))
    return warnings


def process_file(file_path: str, class_name: str, method_name: str) -> bool:
    with open(file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    # Step 1: Find the impl block for the class
    # We look for a line like "impl ClassName" or "impl PyClassName" or "impl<'a> ClassName"
    impl_idx = -1
    for i, line in enumerate(lines):
        if (
            "impl" in line
            and (class_name in line or f"Py{class_name}" in line)
            and "struct" not in line
            and "enum" not in line
        ):
            # Let's verify it's a block opening or has a pymethods attribute near it
            # We can search backwards for #[pymethods] or #[cfg_attr(..., pymethods)]
            is_pymethod = False
            for j in range(max(0, i - 3), i + 1):
                if "pymethods" in lines[j]:
                    is_pymethod = True
                    break
            if is_pymethod:
                impl_idx = i
                break

    if impl_idx == -1:
        return False

    # Step 2: Find the method within the impl block
    # We find "fn method_name" or "pub fn method_name"
    method_idx = -1
    brace_depth = 0
    for i in range(impl_idx + 1, len(lines)):
        line = lines[i]
        # Keep track of braces to stay within the impl block
        brace_depth += line.count("{") - line.count("}")
        if brace_depth < 0:
            # We exited the impl block
            break

        # Check if this line is our method definition
        # Must match fn method_name( or pub fn method_name( or fn method_name < or py_ prefixed equivalents
        if re.search(r"\bfn\s+(py_)?" + re.escape(method_name) + r"\b", line):
            method_idx = i
            break

    if method_idx == -1:
        return False

    # Step 3: Extract return type from the method signature
    # Signature might span multiple lines, so let's combine lines from method_idx forward until we find ')' or '{' or '->'
    sig = ""
    for i in range(method_idx, len(lines)):
        sig += lines[i]
        if "{" in lines[i]:
            break

    # Extract return type
    # e.g., "fn method_name(...) -> ReturnType {"
    # We find the part after "->" and before "{" or "where"
    return_type = ""
    if "->" in sig:
        part = sig.split("->")[1].strip()
        # strip block opening { or where clauses
        if "{" in part:
            part = part.split("{")[0].strip()
        if "where" in part:
            part = part.split("where")[0].strip()
        return_type = part.strip()
    else:
        return_type = "()"

    py_return_type = rust_to_py_type(return_type, class_name)
    print(
        f"Mapped {class_name}.{method_name} (Rust: {return_type}) -> Python: {py_return_type}"
    )

    # Step 4: Find the docstring/decorators above the method and insert :rtype:
    # We search backwards from method_idx to find docstrings starting with "///"
    doc_lines_indices = []
    decorator_lines_indices = []
    for i in range(method_idx - 1, impl_idx, -1):
        line = lines[i]
        if line.strip().startswith("///"):
            doc_lines_indices.append(i)
        elif line.strip().startswith("#["):
            decorator_lines_indices.append(i)
        elif not line.strip():
            # empty line, keep searching
            continue
        else:
            # stop if we hit something else (unless we already found doc lines)
            if doc_lines_indices:
                break
            # if we found only decorators and then something else, there might be no docstring
            if not line.strip().startswith("#[") and not line.strip().startswith("///"):
                break

    if doc_lines_indices:
        # We have a docstring!
        # doc_lines_indices is in reverse order (closest first).
        # Let's find the last (bottom-most) doc line
        last_doc_idx = min(
            doc_lines_indices
        )  # since we searched backwards, min() is the highest line index, which is the last line of consecutive doc comments
        # Wait, actually, let's reverse search logic:
        # Say method is at index 100.
        # Lines 95, 96, 97 are ///. Last doc idx is 97.
        # Let's verify:
        # Indices in doc_lines_indices would be [97, 96, 95] (if consecutive).
        # Bottom-most (last) line of the doc comment is 97, which is max(doc_lines_indices)!
        last_doc_idx = max(doc_lines_indices)

        # Check if :rtype: is already in the docstring (it shouldn't be, but let's check)
        for idx in doc_lines_indices:
            if ":rtype:" in lines[idx]:
                print(f"  :rtype: already present in {class_name}.{method_name}")
                return False

        # Construct indent of the doc line
        indent = lines[last_doc_idx].split("///")[0]
        # Insert "/// :rtype: <py_return_type>" after last_doc_idx
        new_line = f"{indent}/// :rtype: {py_return_type}\n"
        lines.insert(last_doc_idx + 1, new_line)
        print(f"  Inserted rtype at line {last_doc_idx + 2} of {file_path}")
    else:
        # No docstring! Let's insert a docstring right before the method or decorators
        insert_idx = method_idx
        if decorator_lines_indices:
            insert_idx = min(decorator_lines_indices)

        indent = (
            lines[insert_idx].split("#")[0]
            if decorator_lines_indices
            else lines[insert_idx].split("fn")[0]
        )
        new_lines = [f"{indent}/// :rtype: {py_return_type}\n"]
        lines[insert_idx:insert_idx] = new_lines
        print(f"  Created new docstring at line {insert_idx + 1} of {file_path}")

    # Save changes
    with open(file_path, "w", encoding="utf-8") as f:
        f.writelines(lines)
    return True


def find_and_process(class_name: str, method_name: str) -> bool:
    # Search in nyx-core and nyx-py
    search_dirs = ["../nyx-core", "."]
    for search_dir in search_dirs:
        for root, dirs, files in os.walk(search_dir):
            # Ignore some directories
            if (
                "target" in root
                or ".git" in root
                or "node_modules" in root
                or ".venv" in root
            ):
                continue
            for file in files:
                if file.endswith(".rs"):
                    file_path = os.path.join(root, file)
                    if process_file(file_path, class_name, method_name):
                        return True
    return False


def main():
    for module in ["nyx_space.mission_design", "nyx_space.orbit_determination"]:
        print(f"Processing module: {module}")
        warnings = get_missing_rtypes(module)
        print(f"Found {len(warnings)} missing return types in {module}")
        success_count = 0
        for class_name, method_name in warnings:
            if find_and_process(class_name, method_name):
                success_count += 1
            else:
                print(f"Failed to find or annotate {class_name}.{method_name}")
        print(
            f"Successfully annotated {success_count} / {len(warnings)} methods in {module}"
        )


if __name__ == "__main__":
    main()
