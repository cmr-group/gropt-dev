## July 3rd 2026

Easy first steps:
- XX Allow defining the diffusion timings and dt
- XX Output of feasibility status

Necessary steps:
- XX Loading asc files and using those parameters
- XX Multiple SAFE constraints
- XX Allow gropt get_SAFE to take .asc files
- XX Function to fake axis

Helpful:
- Internal TE finder or max b-value finder

Solver Improvements:
- Try NLCG for b-value optimization
- Add simple tik reg into general ILS class (just add lam*identity to lhs)
- TV minimization
    - FISTA/ADMM solver instead of CG
    - Add it to list of constraints, don't do checks


## Next Steps
- Integration with other sequence elements
    - Start with hard-coded blips for crushers and EPI-prewinder
    - Maybe a dummy EPI waveform
    - Test with lower dt for this
    - Full sequence waveforms
- Multi-axis support
    - Test on each gradient axis