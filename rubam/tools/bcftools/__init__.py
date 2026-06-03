"""Pure-Rust ports of bcftools subcommands.

This package wraps the Rust implementations in `src/tools/bcftools/`. Each
function takes the same arguments as the corresponding bcftools subcommand
whenever possible. Implementations are filled in across tasks B2..B8.
"""

from __future__ import annotations

from rubam._rubam import bcftools_view as _view
from rubam._rubam import bcftools_sort as _sort_bcf
from rubam._rubam import bcftools_index as _bcftools_index
from rubam._rubam import bcftools_concat as _concat
from rubam._rubam import bcftools_query as _query
from rubam._rubam import bcftools_stats as _stats
from rubam._rubam import bcftools_norm as _norm


def view(
    input: str,
    *,
    region: str | None = None,
    samples: list[str] | None = None,
    output: str | None = None,
    output_type: str = "v",
    header_only: bool = False,
    no_header: bool = False,
) -> int:
    """bcftools view equivalent. Returns the number of records written."""
    return _view(
        input,
        region=region,
        samples=samples,
        output=output,
        output_type=output_type,
        header_only=header_only,
        no_header=no_header,
    )


def sort(input: str, output: str, *, output_type: str = "v") -> int:
    """bcftools sort equivalent. In-memory sort by (chrom, pos).

    Sorts records by the contig declaration order in the header (not
    alphabetically), then by position. Returns the number of records written.

    output_type: 'v' plain VCF (default), 'z' BGZF VCF.gz, 'b' BCF.
    """
    return _sort_bcf(input, output, output_type=output_type)


def index(input: str, *, csi: bool = False, force: bool = False) -> str:
    """bcftools index equivalent. Returns the path of the written index."""
    return _bcftools_index(input, csi=csi, force=force)


def query(input: str, format: str, *, output: str) -> int:
    """bcftools query equivalent. Returns the number of records emitted.

    The format string supports the standard bcftools placeholders:
      %CHROM, %POS, %ID, %REF, %ALT, %QUAL, %FILTER, %INFO/KEY,
      and per-sample %SAMPLE, %GT, %FORMATKEY inside [ ... ].
    Escape sequences \\t and \\n are recognized.
    """
    return _query(input, format, output)


def concat(inputs: list[str], output: str, *, output_type: str = "v") -> int:
    """bcftools concat equivalent. Glue multiple sorted VCFs / BCFs.

    Headers must be compatible: same sample names and same contig list.
    Returns the total number of records written.
    """
    return _concat(list(inputs), output, output_type=output_type)


def stats(input: str) -> dict:
    """bcftools stats equivalent (subset).

    Returns a dict with:
      total_records, snps, indels, mnps, complex,
      transitions, transversions, ts_tv_ratio,
      samples: { name: {hom_ref, het, hom_alt, missing} }
    """
    return _stats(input)


def norm(
    input: str,
    output: str,
    *,
    output_type: str = "v",
    multiallelic: str = "",
    reference: str | None = None,
) -> dict:
    """bcftools norm equivalent.

    multiallelic='-' splits multi-allelic sites into one ALT per record.
    reference=path enables left-alignment of indels (.fai auto-built).

    Returns a dict: {records_in, records_out, left_aligned}.
    """
    return _norm(
        input,
        output,
        output_type=output_type,
        multiallelic=multiallelic,
        reference=reference,
    )


__all__ = ["concat", "index", "norm", "query", "sort", "stats", "view"]
