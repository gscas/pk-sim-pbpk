# R API Reference for OSPSuite

## Installation

```r
# Install from GitHub
remotes::install_github("Open-Systems-Pharmacology/OSPSuite.R")
# Requires R >= 4.1 and .NET Core runtime
```

## Basic Workflow

```r
library(ospsuite)

# Load simulation from PKML
sim <- loadSimulation("path/to/simulation.pkml")

# Get and set parameter values
lipophilicity <- getParameter("Compound1|Lipophilicity", sim)
setParameterValue(lipophilicity, 2.5)

fu <- getParameter("Compound1|Fraction unbound (plasma, reference value)", sim)
setParameterValue(fu, 0.1)

# Run simulation
results <- runSimulation(sim)

# Get output paths
outputPaths <- getAllOutputPaths(results)
# e.g., "Organism|PeripheralVenousBlood|Compound1|Plasma (Peripheral Venous Blood)"

# Extract time and concentration
time <- results$Time
conc <- getOutputValues(results, outputPath)
```

## Population Simulation

```r
# Create population from PK-Sim
pop <- createPopulation(
  numberOfIndividuals = 100,
  seed = 123,
  populationCharacteristics = PopulationCharacteristics(
    species = Species$Human,
    population = HumanPopulation$European,
    age = ParameterDistribution(mean = 30, deviation = 10, distribution = DistributionTypes$Normal),
    weight = ParameterDistribution(mean = 70, deviation = 15, distribution = DistributionTypes$Normal)
  )
)

# Run population simulation
popResults <- runSimulation(sim, population = pop)

# Extract individual results
for (i in 1:popResults$count) {
  individualResult <- popResults$individualResults[[i]]
  # Process each individual
}
```

## Parameter Identification

```r
# Load observed data
obsData <- loadDataRepositoryFromExcel("observed_data.xlsx")

# Define optimization
optimParams <- c("Compound1|Lipophilicity", "Compound1|Specific clearance")

# Run identification
identResult <- identifyParameters(
  simulation = sim,
  observedData = obsData,
  parametersToIdentify = optimParams,
  optimizationMethod = OptimizationMethod$NelderMead
)
```

## Sensitivity Analysis

```r
saResult <- sensitivityAnalysis(
  simulation = sim,
  outputPath = "Organism|PeripheralVenousBlood|Compound1|Plasma (Peripheral Venous Blood)",
  parametersToVary = c("Compound1|Lipophilicity", "Compound1|Fraction unbound (plasma, reference value)"),
  variationRange = 0.1  # +-10%
)
```

## Export

```r
# Export to PKML
exportSimulationToPKML(sim, "output/simulation.pkml")

# Export results to CSV
exportResultsToCSV(results, "output/results.csv")

# Export snapshot
exportSnapshot(sim, "output/snapshot.json")
```
