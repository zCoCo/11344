from astroquery.mast import Observations
from astropy.timeseries import TimeSeries
from astropy.io import ascii
import os

import csv
from lightkurve import search_lightcurvefile
import time
from os import path

STARTING_INDEX = 1382;  # Used for resuming loads

period_file = './data/transit_periods.csv';
download_cache_dir = "O:/lightkurve_cache";

# Whether to Recollect the Periods for all Stars:
RECOLLECT_PERIODS = True;

with open('./data/star_data.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',');
    count = -1;  # Number of objects loaded (or attempted) already
    start = 0;
    for row in readCSV:
        if count == 0:
            start = time.time();  # don't count setup time in ETA calculation
        if count >   -1 and count >= STARTING_INDEX:  # ignore first row of headers
            tic = row[1];
            if True:
                DT = time.time() - start;
                ETA = (3454-count) * DT / ( (count-STARTING_INDEX) if count > STARTING_INDEX else 1 );
                print("Fetching Light Curve " + str(count) + " for TIC " + str(tic));
                print("Time Elapsed: " + str(int(DT/60)) + "min, ETA: " + str(int(ETA/60)) + "min");
                print(" . . . ");

            lightcurve_file = './light_curves/lc_'+tic+'.csv';
            transit_file = './transits/tr_'+tic+'.csv';
            if RECOLLECT_PERIODS or (not path.exists(lightcurve_file) or not path.exists(transit_file)):
                # Collect all TESS light curves (event observations) for this object:
                if path.exists(download_cache_dir):
                    lcf = search_lightcurvefile('TIC '+tic).download_all(download_dir=download_cache_dir);
                else:
                    lcf = search_lightcurvefile('TIC '+tic).download_all();
                if lcf is not None:
                    try:
                        # Use stitch to combine all observed events into single curve (all points are in absolute time):
                        lcf = lcf.PDCSAP_FLUX.stitch();

                        # Post process data:
                        # Clean up photometric data (remove nans, long term trends, and noise):
                        lcf = lcf.remove_nans().flatten().bin(binsize=5);

                        # Save combined cleaned curve as csv:
                        if not path.exists(lightcurve_file):
                            lcf.to_csv(lightcurve_file);
                            print("Light Curve Saved.");

                        # Fold Data to Find Transit Curve:
                        # Find highest power frequency (.: period)
                        pg = lcf.to_periodogram(oversample_factor=2);
                        period = pg.period_at_max_power;

                        # Save Period to Period File:
                        if RECOLLECT_PERIODS or not path.exists(transit_file):
                            with open(period_file, mode='a') as pf:
                                csv_writer = csv.writer(pf, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL);
                                csv_writer.writerow([str(tic), str(period)]);

                        # Finish Processing Transit Curve:
                        if not path.exists(transit_file):
                            # Fold at period and remove noise:
                            transit = lcf.fold(period=period).bin(binsize=5);
                            #transit = lcf.fold(period=period).bin(binsize=35);

                            # Save Folded Transit Event
                            transit.to_csv(transit_file);
                            print("Transit Curve Saved.");

                    except:
                        print("-----##### Couldn't post-process lightcurve for TIC #####-----" + str(tic));
            else:
                print(" - Skipping - ");
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
