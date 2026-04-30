"""Dev-only script. Emits the markdown table used in the README's
'Mechanisms' section, sourced from `tecap.constants.MECHANISM_DEFINITIONS`.
Run after editing the dict and paste the output into README.md.

Usage:
    python scripts/render_mechanisms.py
"""

from tecap.constants import (
    BUCKET_DEFINITIONS,
    BUCKETS,
    CATEGORIES,
    MECHANISM_DEFINITIONS,
)


def main():
    print("| Bucket | What it means | Why it matters |")
    print("|---|---|---|")
    for cat in CATEGORIES:
        d = MECHANISM_DEFINITIONS[cat]
        print(f"| **{d['short']}** | {d['what']} | {d['why']} |")
    print()
    print("# Basecomp buckets")
    print("| Bucket | What it means | Why it matters |")
    print("|---|---|---|")
    for b in BUCKETS:
        d = BUCKET_DEFINITIONS[b]
        print(f"| **{b}** | {d['what']} | {d['why']} |")


if __name__ == "__main__":
    main()
