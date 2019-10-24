from astroquery.mast import Observations
from astropy.timeseries import TimeSeries
from astropy.io import ascii
import os

import csv
from lightkurve import search_lightcurvefile
import time

with open('./data/star_data.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',');
    count = 0;
    start = time.time();
    for row in readCSV:
        count += 1;
        if count > 1:
            tic = row[1];
            if True:
                # TODO:
                # append all light curves together (point order doesn't matter & they're all in abs time)
                DT = time.time() - start;
                ETA = (3454-count) * DT / ( (count-1) if count > 1 else count );
                print("Fetching Light Curve " + str(count) + " for TIC " + str(tic));
                print("Time Elapsed: " + str(int(DT/60)) + "min, ETA: " + str(int(ETA/60)) + "min");
                print(" . . . ");
            lcf = search_lightcurvefile('TIC '+tic).download();
            if lcf is not None:
                lcf.PDCSAP_FLUX.to_csv('./light_curves/lc_'+tic+'.csv');

'''
Old MAST method:
obs = Observations.query_criteria(target_name="77111228", dataproduct_type="timeseries");
data = Observations.get_product_list(obs);
if len(data) > 0:
    # download data:
    manifest = Observations.download_products(data);
    # find and load data:
    for f in manifest['Local Path']:
        if f.endswith('_lc.fits'): # contains light curve (lc) data
            ts = TimeSeries.read(f, format='tess.fits')
        else: # extraneous file, not needed.
            os.remove(f)
            '''
