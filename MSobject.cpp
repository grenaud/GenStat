/*
 * MSobject
 * Date: Nov-02-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "MSobject.h"

MSobject::MSobject(vector<string> * inputString){ 

    int currentInd=0;
    numberIndividual=int(inputString->size())-2;

    for(int i=0;i<int(inputString->size());i++){
	// cout<<"line "<<inputString->at(i)<<endl;

	if(i==0){   //segsites:
	    if(!strBeginsWith(inputString->at(i),"segsites:")){
		cerr<<"First line should begin with the number of segregating sites line="<<inputString->at(i)<<endl;
		exit(1);
	    }
	    vector<string> temp = allTokens(inputString->at(i),' ');
	    numberSegSites=destringify<int>(temp[1]);
	    if(numberSegSites == 0){
		cerr<<"No segregating sites at line ="<<inputString->at(i)<<endl;
		exit(1);		
	    }
	    alleleMatrix = new bool * [numberIndividual];
	    for(int n=0;n<numberIndividual;n++){
		alleleMatrix[n] = new bool [numberSegSites];
	    }
		
	    
	    continue;
	}

	if(i==1){   //positions:
	    if(!strBeginsWith(inputString->at(i),"positions:")){
		cerr<<"Second line should begin with the number of positions line="<<inputString->at(i)<<endl;
		exit(1);
	    }
	    
	    vector<string> temp = allTokens(inputString->at(i),' ');
	    if( (temp.size()-2) != numberSegSites){
		cerr<<"Wrong number of positions="<<inputString->at(i)<<" expected "<<numberSegSites<<" found  "<<(temp.size()-2) <<endl;
		exit(1);
	    }

	    for(int fieldPos=1;fieldPos<=numberSegSites;fieldPos++){
		//positions.push_back( destringify<double>(temp[fieldPos]) );
		positions.push_back( destringify<unsigned int>(temp[fieldPos]) );
	    }
	    
	    if( (positions.size()) != numberSegSites){
		cerr<<"Wrong number of positions="<<positions.size()<<" expected "<<numberSegSites <<endl;
		exit(1);
	    }

	    continue;	    
	}



	//put bool matrix
	if((inputString->at(i).length()) != numberSegSites){
	    cerr<<"Wrong nunber of site  in "<<inputString->at(i)<<" expected "<<numberSegSites<<" found  "<<(inputString->at(i).size()) <<endl;
	    exit(1);
	}
	
	for(int k=0;k<int( inputString->at(i).length() );k++){
	    if(inputString->at(i)[k] == '0'){
		alleleMatrix[currentInd][k]=false;
	    }else{
		if(inputString->at(i)[k] == '1'){
		    alleleMatrix[currentInd][k]=true;
		}else{
		    cerr<<"Wrong character  in "<<inputString->at(i) <<endl;
		    exit(1);
		}
	    }

	}

	currentInd++;

    }
    
    // for(int n=0;n<numberIndividual;n++){
    // 	for(int m=0;m<numberSegSites;m++)
    // 	    cout<<alleleMatrix[n][m];
    // 	cout<<endl;	    
    // }

}

MSobject::~MSobject(){
    for(int n=0;n<numberIndividual;n++){
      	delete [] alleleMatrix[n];
    }
     delete [] alleleMatrix;
}




const bool * MSobject::getBoolArrayIndividual(int indexInd) const{
    if(indexInd == 0 ){
	cerr<<"Cannot return individual zero, the index is 1-based"<<endl;
	exit(1);
    }

    if(indexInd > numberIndividual ){
	cerr<<"Error: MSobject has "<<numberIndividual<<" cannot return individual #"<<indexInd<<endl;
	exit(1);
    }

    return alleleMatrix[indexInd -1]; //C++  0-based 
}


unsigned int MSobject::getNumberSegSites() const{
    return numberSegSites;
}

int MSobject::getNumberIndividuals() const{
    return numberIndividual;
}

//const vector<double> * MSobject::getPositions() const{
const vector<int> * MSobject::getPositions() const{
    return &positions;
}
