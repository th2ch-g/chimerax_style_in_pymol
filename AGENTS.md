# Development

- Reply in concise Japanese. Use English for code, documentation, and commits.
- Run Python through uv with an interpreter that provides PyMOL.
- Keep all paths portable. Keep runtime independent of ChimeraX and sibling repositories.
- Read docs/pymol_chimerax.md before changing the public interface.
- Preserve source coordinates, atom properties, and unrelated representations.
- Validate geometry in real PyMOL, including GPU and ray images in a separate process.
- Keep downloaded inputs, environments, caches, and temporary renders ignored.
- Published docs/gallery PNGs are the exception: version them for README display.
- Keep the English and Japanese guides synchronized.
- Scope searches to relevant paths. Retain attribution in NOTICE.
- Append a blank line and Co-authored-by: Codex <noreply@openai.com> to commits.
