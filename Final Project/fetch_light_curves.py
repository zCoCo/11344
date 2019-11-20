from astroquery.mast import Observations
from astropy.timeseries import TimeSeries
from astropy.io import ascii
import os

import csv
from lightkurve import search_lightcurvefile
import time

STARTING_INDEX = 107;  # Used for resuming loads

with open('./data/star_data.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',');
    count = -1;  # Number of objects loaded (or attempted) already
    start = 0;
    for row in readCSV:
        if count == 0:
            start = time.time();  # don't count setup time in ETA calculation
        if count > -1 and count >= STARTING_INDEX:  # ignore first row of headers
            tic = row[1];
            if True:
                DT = time.time() - start;
                ETA = (3454-count-STARTING_INDEX) * DT / ( count-STARTING_INDEX if count > STARTING_INDEX else 1 );
                print("Fetching Light Curve " + str(count) + " for TIC " + str(tic));
                print("Time Elapsed: " + str(int(DT/60)) + "min, ETA: " + str(int(ETA/60)) + "min");
                print(" . . . ");
            # Collect all TESS light curves (event observations) for this object:
            lcf = search_lightcurvefile('TIC '+tic).download_all();
            if lcf is not None:
                try:
                    # Use stitch to combine all observed events into single curve (all points are in absolute time):
                    lcf = lcf.PDCSAP_FLUX.stitch();

                    # Post process data:
                    # Clean up photometric data (remove nans, long term trends, and noise):
                    lcf = lcf.remove_nans().flatten().bin(binsize=5);

                    # Save combined cleaned curve as csv:
                    lcf.to_csv('./light_curves/lc_'+tic+'.csv');

                    # Fold Data to Find Transit Curve:
                    # Find highest power frequency (.: period)
                    pg = lcf.to_periodogram(oversample_factor=1);
                    period = pg.period_at_max_power;
                    # Fold at period and remove noise:
                    transit = lcf.fold(period=period).bin(binsize=5);

                    # Save Folded Transit Event:
                    transit.to_csv('./transits/tr_'+tic+'.csv');

                except:
                    print("-----##### Couldn't post-process lightcurve for TIC #####-----" + str(tic));
        count += 1;

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
