/*
 *  CBDataFromPressureFile.cpp
 *  CardioMechanics
 *
 *  Created by ek717 on 25.04.2020
 *  Copyright 2011 Institute for Biomedical Engineering, Karlsruhe Institute of Technology (KIT). All rights reserved.
 *
 */

#include "CBDataFromPressureFile.h"
#include <algorithm> // needed for std::min()
#include <stdlib.h> // needed for realpath()
#include "filesystem.h"
#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>



CBDataFromPressureFile::CBDataFromPressureFile()
{
    dataTimeSteps_.clear();
    dataPressureAtTimeStep_.clear();
}


CBDataFromPressureFile::CBDataFromPressureFile(ParameterMap* parameters, std::string parameterKey)
{
    CBDataFromPressureFile();
    
    std::string filename = parameters->Get<std::string>(parameterKey + ".Filename");
    
    Init(filename);
}


void CBDataFromPressureFile::Init(std::string filename)
{
    // LoadFileList(filename);
    LoadDataSet(filename); //0);
}





bool CBDataFromPressureFile::LoadDataSet(std::string filename)
{

    int lineCtr = 0;
    bool reachedEOF = false;
    std::vector<std::string> temp;
    
    std::string line;
    std::ifstream ifs;
    ifs.open(filename.c_str(), std::ifstream::in);    //open file to read
    
    //check if file is ok
    if(!ifs.good()){
        
        throw std::runtime_error("CBDataFromPressureFile::LoadDataSet(std::string filename): File list: " + filename + " does not exist !!!");
        ifs.close();
        return 0;
    } else {
        
        while (ifs.good()) {
            getline(ifs, line);    //read next line from file, save in: line
            lineCtr++;
            
            if (ifs.eof() && line.empty()) {
                //reached end of file, which is a empty line. Return with success code 1
                ifs.close();
                return 1;
            } else if (ifs.eof() && !line.empty()) {
                //reached end of file in last line of document. Set reachedEOF to true,
                //so this last line gets processed and after that close ifs
                reachedEOF = true;
            }
            
            //if (lineCtr == 1) {
                //the first line of the file contains header information. So skip this line or
                //maybe use the information for error prevention (e.g. too few/much data points given, ...)
               // continue;
            //}
            
            //process new dataset
            temp = SplitString(line);
            
            if (temp.size() != 2) {
                //if it doesn't hold 2 elements, something went wrong. Return with error code 0.
                ifs.close();
                return 0;
            }
   
            // add data
            dataTimeSteps_.push_back(stof(temp.at(0)));
            dataPressureAtTimeStep_.push_back(stof(temp.at(1)));
            
            //check whether it was the last line and it contained also the EOF
            if (reachedEOF) {
                //EOF and last data item were in the same line. So close ifs and return 1 for success.
                ifs.close();
                return 1;
            }
        }//end while
    }//end if

    
    return(1);
    
}


TFloat CBDataFromPressureFile::Get(TFloat time, TInt index)
{
    
    TFloat pressureCurrent = 0;
    bool timeInPressureDataFile = false;
    TInt ind = 0;
    
    while(time > dataTimeSteps_[ind] )
    {
        ind++;
        timeInPressureDataFile = true;
        // time is bigger then the last time step from the input pressure file
        if(ind > dataTimeSteps_.size()){
            timeInPressureDataFile = false;
            ind--;
            break;
        }
    }

    if(timeInPressureDataFile){
        if(ind == dataTimeSteps_.size()){
            if(dataTimeSteps_.size() == 1){
                // pressure file contains only one time step
                pressureCurrent = dataPressureAtTimeStep_[ind];
                return pressureCurrent;
            }/*else{
                //  to interpolate between one before last and last time step shift the index to the step before
                ind--;
                    
            }*/
        }
        ind--;
        if(dataTimeSteps_[ind] == dataTimeSteps_[ind+1])
        throw std::runtime_error("CBDataFromPressureFile::Get(std::string filename): In the pressure input file there is the same timestep twice !!!");
        // interpolate between current and next time step
        pressureCurrent = (dataPressureAtTimeStep_[ind] * (dataTimeSteps_[ind+1] - time) + dataPressureAtTimeStep_[ind+1] *  (time - dataTimeSteps_[ind])) / (dataTimeSteps_[ind+1]-dataTimeSteps_[ind]) ;
            
        
    }
        
    
    
    return pressureCurrent;
}



/// Method to split an string at whitespaces.
std::vector<std::string> CBDataFromPressureFile::SplitString(std::string line) {
    
    std::vector<std::string> tokens;
    
    char delimiter = ' ';
    std::string tok;
    std::string::iterator it=line.begin();
    while (  it!=line.end() ) {
        if (*it != delimiter && *it != '\r'){
            tok.push_back(*it);
        } else {
            if (!tok.empty()){
                tokens.push_back(tok);
                tok.clear();
            }

        }
        ++it;
    }
    // if there is a carriage return ('\r') as a newline character, tok is empty and we don't need to push it into token.
    if (tok.compare("") != 0) {
        tokens.push_back(tok);
        tok.clear();
    }
    // lh326: Note there should be an empty line at the end of the pressureGradients.txt (p.txt)!!
    return tokens;
}

