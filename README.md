# xrs_example
Some examples of working with GOES/XRS data:

- To get data look at https://github.com/ianan/xrs_example/blob/main/xrs15_example.ipynb
- To calculate the T,EM look at https://github.com/ianan/xrs_example/blob/main/xrs15_tem_example.ipynb

#### Note:
- Best to work with the avg1min data for all GOES sat as removes unwanted things like electron contamination, particularly if looking at A/B Class levels. Though this data is only 1min cadence and isn't at the moment available via sunpy/fido: can get the daily files directly from NOAA, i.e. for GOES16 https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/xrsf-l2-avg1m_science/
