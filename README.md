# UVa decoder
A modified version of the UVa MPD4 decoder. 

## main change

* changed the bank number and the tag 
* change the data search structure of the 'inputhandler.cpp' so as to decode the PRex data 
* changed the grephic out style

## Dev Plan
### Raw display function 

- [ ] Change the Raw Display
- [ ] in the raw display add the cluster searching 
- [ ] add the text information on the canvas 
- [ ] display the frame size and frame check information


## pedestal generator instructions
### 1. check the whether the mapping match the prex mapping
### 2. change the config/gem.cfg file 
    # runType
    RUNTYPE: CalPedestal
    # pedestal
    SAVEPED: ./Pedestal/     // this is the file name of the generated pedestal
    # # Input File for physics analysis; NOT FOR OTHER TYPE ANALYSIS

    INPUTFILE:   //  the pedestal raw run file name


### 3. generate the pedestal
    ./mpd4_decoder

### 4. save the file to the prex database
* it will generate a file named PRex_Pedestal.txt
* this file is conpatible with the standard prex pedestal, need to replace the pedestal in the prex database


## RUN LIST

#### 16-Jun-2019
* Pedestal Run 1371 20513

#### 23-Aug-2019
* Pedestal Run 2155, 21297
   
#### 03-Sep-2019
* Pedestal Run 2303, 21425
* HV off
* pulse trigger
