# Private sample production on CRAB

Follow the installation instructions of the main [README.md](../README.md) 

Make sure to activate the proxy
```
voms-proxy-init --voms cms --valid 186:00
```
Then, 
```
cd HNLsGen/crab
```
List the desired mass/ctau points to generate in `<points>.py`

Make sure to have the necessary fragments/pdl files copied in the respective directories.

Run the tool as
```
python CRABLauncher.py --pl <prod_label> --points <points>.py --nevents <nevents> <--dosubmit>
```
with nevents the requested number of events at analysis level.

Note that the number of jobs is adapted such that each job processes 100k events.
