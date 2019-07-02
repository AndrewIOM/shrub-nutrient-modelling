Shrub Nutrient Modelling
======

Mechanistic models of Arctic shrub growth, accounting for nitrogen and temperature limitation. 

This library contains [Bristlecone](https://github.com/AndrewIOM/bristlecone) scripts that have been used in three upcoming papers.

Requirements:
- .NET Core (2.1 or higher)
- Paket installed as a global tool (`dotnet install tool Paket -g`)
- For R ggplot traceplots, you must currently run F# interactive using mono or the full .NET framework, with a version of R below 3.50.

Steps to use:
1. Download Paket packages: `dotnet paket install`.
2. Run scripts line-by-line interactively using F# interactive: `dotnet fsi`.