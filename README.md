# Power-Analysis-for-PAPE-Effects-in-Bench-Press-vs-Military-Press-on-Isometric-Performance

📌 Project Overview

This project uses simulation-based power analysis to estimate whether a typical sports science study design is sufficiently powered to detect Post-Activation Performance Enhancement (PAPE) effects in upper-body performance.

The analysis compares bench press (BP) and seated military press (SMP) as conditioning activities and their impact on isometric strength performance (IBP, ISMP).

❓ Research Question

Is a sample size of N = 30 sufficient to reliably detect a small-to-moderate PAPE effect (Cohen’s d ≈ 0.3) in isometric upper-body performance following BP and SMP conditioning?

📚 Background

Post-Activation Performance Enhancement (PAPE) refers to acute improvements in muscular performance following high-intensity voluntary contractions.

Previous literature suggests:

Typical effect sizes around d ≈ 0.15–0.5 (commonly ≈ 0.3)
Most studies focus on bench press and ballistic outcomes
Limited evidence exists for:
Isometric performance outcomes
Vertical pushing movements (e.g., military press)
🧪 Study Design (Simulated)
Planned sample size: N = 30 participants
Within-subject design:
Control vs PAPE condition
Two conditioning activities: BP and SMP
Multiple time points: 0, 5, 7, 10 minutes
Measured outcomes:
Isometric Bench Press (IBP)
Isometric Seated Military Press (ISMP)
📊 Simulation Approach
Monte Carlo simulation (10,000 iterations)
Effect size assumption: d ≈ 0.3

Data generated based on:

Pilot measurements
Literature-based estimates

Statistical model:
Linear mixed-effects models (LMM)

Fixed effects:

condition (control vs PAPE)
time point
interaction (condition × time)

Random effects:

participant ID
📈 Key Results

Estimated statistical power for detecting PAPE effect (d ≈ 0.3):

Power ≈ 0.99 (99%) with N = 30

Additional analysis:

Equivalence testing power ≈ 0.78 (78%) for comparing IBP vs ISMP effects

Power estimation based on:

proportion of simulations where 95% CI excluded zero
🧠 Interpretation

The simulation suggests that:

A sample size of N = 30 is sufficient to detect small-to-moderate PAPE effects in isometric strength outcomes
Bench press and military press conditioning are expected to produce comparable PAPE effects
Within-subject designs enable high statistical power even with moderate sample sizes

This has important implications for:

evaluating movement-specificity of PAPE
extending PAPE research beyond ballistic tasks
improving experimental design and reproducibility in strength research
🧰 Tools

R, Monte Carlo simulation, linear mixed-effects models (lme4), TOSTER

See the preregistration of the project: https://osf.io/xv9kt/overview
