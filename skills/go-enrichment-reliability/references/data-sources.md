# Data Sources and Tool References

As of 2026-02-25, use these official sources for capability checks and maintenance signals.

## Primary references

- GSEApy GitHub releases: <https://github.com/zqfang/GSEApy/releases>
- clusterProfiler Bioconductor page: <https://bioconductor.org/packages/clusterProfiler/>
- Enrichr API help: <https://maayanlab.cloud/Enrichr/help#api>
- g:Profiler API docs: <https://biit.cs.ut.ee/gprofiler/page/apis>
- WebGestalt API docs: <https://www.webgestalt.org/WebGestalt_2024_Manual_API.html>
- Reactome data release archive: <https://reactome.org/tag/release>
- DAVID update notes: <https://david.ncifcrf.gov/content.jsp?file=updates.html>
- Broad GSEA user guide: <https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm>

## Licensing and access notes

- Enrichr/g:Profiler/WebGestalt/Reactome are web services and can have API limits or uptime variability.
- KEGG has separate licensing terms for some commercial uses; verify usage policy before pathway-scale automation.
- Broad GSEA baseline requires local CLI/runtime setup when strict offline reproducibility is required.

## Operational policy

- Treat web baseline failures as warnings in normal runs.
- Treat web baseline failures as hard failures in `--benchmark-mode`.
- Keep source URLs in reports for traceability.
