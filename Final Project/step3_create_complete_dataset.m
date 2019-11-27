% After fetch_light_curves.py has collected TESS photometry, this will
% call all necessary processes to process the data, add metaparameters, 
% clean the data, post process the results, then partition the dataset.
function step3_create_complete_dataset()
    disp(" #### ADDING TRANSIT FEATURES TO DATASET #### ");
    create_transit_dataset(); % Add transit data to star dataset
    
    disp(" #### ADDING METAPARAMETERS TO DATASET #### ");
    add_metaparameters(); % Add metaparameters obtained from error analysis
    
    disp(" #### POSTPROCESSING DATASET #### ");
    postprocess_data(); % Clean and normalize the data
    
    disp(" #### PARTITIONING DATASET #### ");
    divide_dataset(); % Split Dataset into Training Set, Testing Set, etc.
end