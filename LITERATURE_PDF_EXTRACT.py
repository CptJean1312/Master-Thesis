#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

from pypdf import PdfReader


def iter_pdfs(path: Path) -> Iterable[Path]:
    if path.is_file() and path.suffix.lower() == ".pdf":
        yield path
        return

    if path.is_dir():
        for pdf in sorted(path.rglob("*.pdf")):
            if pdf.is_file():
                yield pdf
        return

    raise FileNotFoundError(f"No PDF file or directory found at: {path}")


def extract_pdf_text(pdf_path: Path) -> str:
    reader = PdfReader(str(pdf_path))
    pages = []
    for idx, page in enumerate(reader.pages, start=1):
        text = page.extract_text() or ""
        pages.append(f"[[PAGE {idx}]]\n{text}")
    return "\n\n".join(pages)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract text from one PDF or all PDFs in a directory."
    )
    parser.add_argument("input_path", help="PDF file or directory containing PDFs")
    parser.add_argument(
        "--output-dir",
        help="Directory for extracted .txt files. If omitted for a single PDF, text is printed to stdout.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional limit on how many PDFs to process from a directory.",
    )
    parser.add_argument(
        "--error-log",
        help="Optional path for a log file listing PDFs that could not be extracted.",
    )
    args = parser.parse_args()

    input_path = Path(args.input_path).expanduser().resolve()
    pdfs = list(iter_pdfs(input_path))

    if args.limit is not None:
        pdfs = pdfs[: args.limit]

    if args.output_dir:
        output_dir = Path(args.output_dir).expanduser().resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        errors = []

        for pdf_path in pdfs:
            try:
                text = extract_pdf_text(pdf_path)
                out_path = output_dir / f"{pdf_path.stem}.txt"
                out_path.write_text(text, encoding="utf-8")
                print(f"Wrote {out_path}")
            except Exception as exc:
                message = f"FAILED\t{pdf_path}\t{type(exc).__name__}: {exc}"
                errors.append(message)
                print(message)

        if args.error_log:
            error_log = Path(args.error_log).expanduser().resolve()
            error_log.parent.mkdir(parents=True, exist_ok=True)
            error_log.write_text("\n".join(errors), encoding="utf-8")

        if errors:
            raise SystemExit(1)
        return

    if len(pdfs) != 1:
        raise SystemExit("Without --output-dir, input_path must point to exactly one PDF.")

    print(extract_pdf_text(pdfs[0]))


if __name__ == "__main__":
    main()
