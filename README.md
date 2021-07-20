# volcanic_event_source
This code maps the spatio-temporal relationships between vents to model eruptive event sources

The development and context of this code are discussed in Gallant et al., (submitted) - this paper is currently in review at Volcanica.

This is a Matlab code developed using MATLAB R2020a. We provide a .xls file (vent_data_all.xls) that contains the data described in the above publication, as well as the outputs described. 

Running the code: 

You need to download this entire directory - all those extra .m files are functions that are called by the code at various points - you need only to download them and stick them in the exact same directory as the 'run_event_source_model.m' file.

There are 4 lines of information you need to customize to run the code - they are described at the top of the 'run_event_source_model.m' file. We provide the required hillshades used to plot the cases described in Gallant et al., 202X in a sub folder in this repository. 

In order to run the code, you simply type 'run_event_source_model' in the matlab command window, or press the triangular green button in the editor tab to run. 

Outputs:

Two figures will be automatically generated, a .fig file of the vents and events and a .fig file of spatial density contours. In order to track down the additional data about the events (the average age of the event, the vents that comprise it, etc -- we show this in the vent_data_all.xls under the ...event coordinates tabs), you need to type in Array = write_event_csv(Events,'####.csv'); --> the #### should be the name of your volcanic field. An automatic list of the event source coordinates is generated as a .csv file --> this csv can be opened in matlab but I have had issues opening in LibreOffice

The sub-folder "additional command function" might help you get some of the outputs you're interested in, so we've included those just in case there is something in there you're interested in. 

