Instructions for use:
  1. Start RStudio.
  
  
  2. Install the required packages:
     a. Locate the "Console" pane in RStudio. It's usually in the bottom left corner.
     b. In that pane, at the prompt that looks like ">", type or paste this command - 
        install.packages(c("readxl", "xlsx", "dplyr", "ggplot2", "lazyeval", "shiny", "shinyjs", "DT"))
     c. Press Enter. This will take a few seconds. You should see several messages telling you that a certain package has been successfully unpacked etc. The last message should say something like -
       The downloaded binary packages are in
	        C:\...\downloaded_packages.
	     Here "..." is a directory path.
	     
	     
	3. Set the default directory:
	   a. Go to the Tools menu.
	   b. Select Global Options.
	   c. Set the "Global working directory (when not in a project):" to Wherever this code has been unzipped, i.e. the folder that has this file.
	   d. Close and then restart RStudio.
	   
Steps 2 and 3 need to be executed only once.

	4. If at least one of "ui.R" or "server.R" from this folder isn't Open in RStudio:
	   a. Go to the File menu.
	   b. Select "Open File...".
	   c. Choose one of those two files from this folder.
	   

	5. Run the app:
	   a. With one of "ui.R" or "server.R" open in RStudio, click on "Run App". It should be a button towards the top right corner of the tab that is showing the file.
	   b. This should pop up a window that shows the user interface(UI).
	   c. To access all features of the UI, click the button at the top of that window, that says "Open in Browser".
	   d. This should open the same UI in your default web-browser e.g Firefox, Chrome, Edge etc.
	   
	   
	6. Use the app:
	   a. Upload the Excel (.xslx) file containing legionella data using the "Browse..." button under "Choose legionella data file".
	   
	   b. While the file is being uploaded and the data is being loaded into the app, you'll see a progress bar and various messages in the bottom right corner of the browser window.
	   
	   c. Until the file has been uploaded, the message "Please choose a legionella data file" will be shown in the tabs on the right and the "Download" buttons will stay disabled. They will stay disabled until the tabs are showing something in their tables.
	   
	   d. After the data has been loaded from the file, "Generate tables" in the bottom left of the window can be clicked to generate the output which will be loaded in the tabs. Before the output is generated, a progress bar and various messages in the bottom right corner of the window will be shown. The data in the tab "Gene fields" is derived from only the specimens, genes and sequences (nodes) shown in the tab "All fields".
	   
	   e. Fields under "Filter By" allow the output to be limited to certain values.
	     I. There are drop-downs for species, serogroups, sequence types and genes. These are all loaded from the uploaded Excel file. All of these allow multiple selections.
	    II. The file "args.xlsx" which is in the same folder as this one, has a column called "Genes". If you'd like the drop-down for genes to have some selected by default, you can specify them in this column. They must all be spelled correctly, separated by commas and have no spaces.
	   III. The fields "Seq Total % Coverage" and "Seq Weighted % Identity" are used to filter the columns with the same name on the tab "Gene fields". These columns show the calculated total % coverage and weighted % identity for each sequence (node). These fields must have valid numbers entered in them (e.g. 20, 20.0, 20.01 etc.). If they don't, the error stating that will be shown under each tab.
	   
	   f. The tables in the tabs have several functionalites:
	     I. The data shown in a table can be further filtered by using the "Search" field in the top right corner of the tab. Any data shown in the table can be entered in this field to limit the display to only that data. 
	    II. The tables are "paged", which means that only a certain range of rows is visible in the table. To view the next or previous range, the buttons at the bottom of the table can be used.
	   III. The range can be altered by the using the drop-down between "Show" and "entries" in the top left corner of the tab.
	   
	   g. The data showing in each tab can be downloaded in the form of an Excel (.xslx) file by clicking the "Download" button in the tab. There are a few things to note here:
	     I. For this to work, the app must be running in a browser as mentioned in step 5. 
	    II. Some columns in the downloaded files don't have the same names as the ones in the tables, because of issues with how R handles the download. So you may need to rename these columns if you want them to match the tables. 
	   III. Only the filtering done using fields mentioned in step 6.e will affect the downloaded data, not the filtering done by "Search" fields mentioned in step 6.f. The latter is only for display purposes.
	  
	  
	  7. Stopping the app:
	    a. Close the browser window or tab showing the app.
	    
	    b. Go back to RStudio.
	    
	    c. Locate a hexagonal "STOP" button right above the "Console" pane mentioned in step 2.a and click it. 
	   
	  8. If the app is running in a browser and the browser's refresh button is clicked, the app may not work correctly. So you will need to restart it. To do so, follow step 8 and then step 5.
	    
	  
	  