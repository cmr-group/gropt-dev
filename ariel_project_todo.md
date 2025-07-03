## July 3rd 2026

Easy first steps:
- XX Allow defining the diffusion timings and dt
- XX Output of feasibility status

Necessary steps:
- Loading asc files and using those parameters
- Multiple SAFE constraints
- Allow gropt get_SAFE to take .asc files
- Function to fake axis

Helpful:
- Internal TE finder or max b-value finder


## Next Steps
- Integration with other sequence elements
    - Start with hard-coded blips for crushers and EPI-prewinder
    - Maybe a dummy EPI waveform
    - Test with lower dt for this
    - Full sequence waveforms
- Multi-axis support
    - Test on each gradient axis