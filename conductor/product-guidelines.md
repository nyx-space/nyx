# Product Guidelines

## Prose Style
All documentation, including API comments, READMEs, and technical guides, should adhere to a **Formal & Technical** style. The tone must be professional, authoritative, and precise, focusing on technical clarity and the accurate use of astrodynamics terminology.

## Branding
The project should be presented with a **Nyx Core Focus**. The primary name is "Nyx," referring to the core astrodynamics engine. When referring to the broader ecosystem or the development organization, use "Nyx Space."

## Documentation Priorities
- **Example-Driven Documentation:** Maintain a comprehensive set of runnable examples for common mission scenarios (e.g., orbit propagation, stationkeeping, orbit determination). These examples serve as both tests and primary learning material for users.
- **Theoretical & Mathematical Guides:** Provide deep-dive documentation that explains the underlying physics, mathematical models, and algorithms implemented within Nyx. This ensures transparency and helps users understand the "why" behind the code.

## API UX Principles
- **Informative Error Handling:** Prioritize clear, actionable, and context-aware error messages. Every error should help the user understand what went wrong and how to fix it, especially in complex mission scenarios.
- **API Discoverability:** Design the library structure and public interfaces to be intuitive and easy to navigate. Use clear naming conventions and logical module hierarchies so users can easily find the tools they need.
- **Consistency:** Ensure that common patterns (e.g., coordinate systems, time formats, data structures) are applied consistently across all modules.