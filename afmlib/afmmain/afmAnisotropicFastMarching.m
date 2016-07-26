function fakeoutput = afmAnisotropicFastMarching(Tvalue, Ttag, F, Boundary, dx, dy, dz, afmSize, D, timeLimit)
% trial list holds the right trial values.
% trialC list holds the values calculated from the group velocity stuff.
% chKnown list holds the known points whose value have changed during iterations.
keySet =   [1000000];
valueSet = [327.2];
trial = containers.Map(keySet,valueSet)
trialC = containers.Map(keySet,valueSet)
N = 0;
fprintf('afmAnisotropicFastMarching function has been called\n');
%The following code creates multidimensional cell of value_iter
%MInitialize matrix for knownChangedImage  
	for(k=1:afmSize(3))
	    % value_iter[k] = new multimap<float,long>::iterator**[Size[1]];
	    % knownChangedImage[k] = new bool*[Size[1]];
	    for(j=1:afmSize(2))
			% value_iter[k][j] = new multimap<float,long>::iterator*[Size[0]];
			% knownChangedImage[k][j] = new bool[Size[0]];
		 	for(i=1:afmSize(1))
		    	% value_iter[k][j][i] = new multimap<float,long>::iterator[2]; % first - trialValue, second - trialCValue
		    	% knownChangedImage[k][j][i] = false;
		    	if(Boundary(k,j,i) == 0)
					N = N + 1;
				end
		    end
		end
	end
	value_iter{afmSize(3), afmSize(2), afmSize(1), 2} = [];
	knownChangedImage = ones(afmSize);
	knownChangedImage = permute(knownChangedImage, [3 2 1]);
	% creating the tetrahedra points.
	% SetTetrahedra(trX, trY, trZ);
	[trX, trY, trZ] =  afmSetTetrahedra();
	% yon[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
	yon = [1,0,0;-1,0,0;0,1,0;0,-1,0;0,0,1;0,0,-1];
	% initial tagging for everypoint in the domain.
  	for k=1:afmSize(3)
    	for j=1:afmSize(2)
	    	for i=1:afmSize(1) 
	      		if(Boundary(k,j,i) == 0 && Ttag(k,j,i) == 200) %  if it's outside the boundary and is known.
		    		for(yonNo=1:6)
				    	x = (i - yon(yonNo,1));
				    	y = (j - yon(yonNo,2));
				    	z = (k - yon(yonNo,3));
				    	conditionone = x>=0 && x<afmSize(1) && y>=0 && y<afmSize(2) && z>=0 && z<afmSize(3);
				    	conditiontwo = Boundary(z,y,x) == 0;
				    	conditionthree = Ttag(z,y,x) ~= 200;
				    	conditionfinal = conditionone && conditiontwo && conditionthree;   
				    	if(conditionfinal)
							fprintf('afmSize(1) is : %d afmSize(2) is : %d afmSize(3) is %d\n', afmSize(1), afmSize(2), afmSize(3));		  
							[Tvalue, Ttag, value_iter] = afmUpdateNeighborhoodTrial(Tvalue, Ttag, F, Boundary, dx, dy, dz, afmSize, D, x, y, z, trial, trialC, trX, trY, trZ, 8, value_iter);
						end
		    		end
				end
	    	end
		end
    end
fakeoutput = []
end

% #include "mxAnisoDistanceTransform.h"
% #include "mex.h"
% #include<math.h>
% using namespace std;
% void AnisotropicFastMarching( float*** Tvalue, int*** Ttag, 
% 							       float*** F, unsigned char*** Boundary, 
% 							       float dx, float dy, float dz, 
% 							       int* Size, float**** D, float timeLimit)
% {
% // trial list holds the right trial values.
% // trialC list holds the values calculated from the group velocity stuff.
% // chKnown list holds the known points whose value have changed during iterations.
%   harita *fm_map;
%   fm_map = new harita();
  
%   //   
%   int trNo,x,y,z,yonNo,i,j,k;
%   int updatedDirection[3];
%   float *p_value;
%   p_value = new float;
%   double *posx,*posy,*posz;
%   posx = new double;
%   posy = new double;
%   posz = new double;
%   long index;
%   time_t  startTrial, endTrial, startKnown, endKnown;
%   double t_find, t_updateTrial, t_updateKnown;
%   /////// TRYING TO GET RID OF FIND STUFF ///////////////////////
%   multimap<float,long>::iterator ****value_iter;
%   bool ***knownChangedImage; // image that keeps the information if a known point changed or not.
  
%   value_iter = new multimap<float,long>::iterator***[Size[2]];
%   knownChangedImage = new bool**[Size[2]];
  
  
%   int N = 0;
%   mexPrintf("size(2): %d, size(1): %d, size(0): %d.\n",Size[2], Size[1], Size[0]);
%   for(k=0;k<Size[2];k++)
%     {
%       value_iter[k] = new multimap<float,long>::iterator**[Size[1]];
%       knownChangedImage[k] = new bool*[Size[1]];
      
%       for(j=0;j<Size[1];j++)
% 	{
% 	  value_iter[k][j] = new multimap<float,long>::iterator*[Size[0]];
% 	  knownChangedImage[k][j] = new bool[Size[0]];
	  
% 	  for(i=0;i<Size[0];i++)
% 	    {
% 	      value_iter[k][j][i] = new multimap<float,long>::iterator[2]; // first - trialValue, second - trialCValue
% 	      knownChangedImage[k][j][i] = false;
% 	      if( Boundary[k][j][i] == 0 )
% 		{
		  
% 		  N++;
% 		}
% 	    }
% 	}
%     }
%   /////////////////////////////////////////////////////////////
  
  
%   float **trX;
%   float **trY;
%   float **trZ;
%   trX = new float*[8];
%   trY = new float*[8];
%   trZ = new float*[8];
%   for(trNo = 0; trNo<8; trNo++)
%     {
%       trX[trNo] = new float[3];
%       trY[trNo] = new float[3];
%       trZ[trNo] = new float[3];
%     }
%   // creating the tetrahedra points.
%   SetTetrahedra(trX, trY, trZ);
  
%   int yon[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
  
  
%   // initial tagging for everypoint in the domain.
%   for(k=0;k<Size[2];k++)
%     {
%       for(j=0;j<Size[1];j++)
% 	{
% 	  for(i=0;i<Size[0];i++)
% 	    {
% 	      if(Boundary[k][j][i] == 0 && Ttag[k][j][i] == 200) //  if it's outside the boundary and is known.
% 		{
% 		  for(yonNo=0;yonNo<6;yonNo++)
% 		    {
% 		      x = (int)(i - yon[yonNo][0]);
% 		      y = (int)(j - yon[yonNo][1]);
% 		      z = (int)(k - yon[yonNo][2]);
% 		      if( x>=0 && x<Size[0] && y>=0 && y<Size[1] 
% 			  && z>=0 && z<Size[2] && Boundary[z][y][x] == 0 
% 			  && Ttag[z][y][x] != 200 )
% 			{
% 			  mexPrintf("first dimension is : %d second dimension is : %d third dimension is %d\n", Size[0], Size[1], Size[2]);		  
% 			  UpdateNeighborhoodTrial( Tvalue, Ttag, F, Boundary, dx, dy, dz, Size, D, x, y, z, fm_map, trX, trY, trZ, 8, value_iter );
			  
% 			}
% 		    }
% 		}
% 	    }
% 	}
%     }
  

%   t_find = 0;
%   t_updateTrial = 0;
%   t_updateKnown = 0;

%   //NOT TO MAKE THE UNNECESSARY SIMULATIONS
%   bool limitReached=false;
  
%   while((fm_map->trial.size() != 0 || fm_map->trialC.size() != 0 
% 	|| fm_map->chKnownX.size() != 0) && !limitReached)
%     {
      
%       if(fm_map->chKnownX.size() == 0)
% 	{
% 	  if(fm_map->trial.size() != 0) // make the minimum element known.
% 	    {
% 	      // finding the minimum element in the good list.
% 	      index = (fm_map->trial.begin())->second;
% 	      ind2sub( posx, posy, posz, (double)index, Size );
% 	      i = lround(*posx);
% 	      j = lround(*posy);
% 	      k = lround(*posz);
% 	      // make the point known and remove it from the good list.
% 	      Ttag[k][j][i] = 200;
% 	      fm_map->trial.erase( fm_map->trial.begin() );
	      
% 	      // NOT TO MAKE THE UNNECESSARY SIMULATIONS.
% 	      if(Tvalue[k][j][i] >= timeLimit)
% 		limitReached=true;
	      
% 	    }
% 	  else // the good list is empty. make the minimum element known. 
% 	    {
% 	      // find the minimum element in the bad list.
% 	      index = (fm_map->trialC.begin())->second;
% 	      ind2sub( posx, posy, posz, (double)index, Size );
% 	      i = lround(*posx);
% 	      j = lround(*posy);
% 	      k = lround(*posz);
% 	      // make the point known and remove it from the bad list.
% 	      Ttag[k][j][i] = 200;
% 	      fm_map->trialC.erase( fm_map->trialC.begin() );
% 	    }
	  
% 	}
%       else // recursive correction of known points.
% 	{
% 	  i = *(fm_map->chKnownX.begin());
% 	  j = *(fm_map->chKnownY.begin());
% 	  k = *(fm_map->chKnownZ.begin());
% 	  fm_map->chKnownX.erase(fm_map->chKnownX.begin());
% 	  fm_map->chKnownY.erase(fm_map->chKnownY.begin());
% 	  fm_map->chKnownZ.erase(fm_map->chKnownZ.begin());
% 	  knownChangedImage[k][j][i] = false; // removed
	  
% 	}
      
      
      
%       // RECURSIVE CORRECTION OF THE KNOWN POINTS //
%       time(&startKnown);
      

      
      

%       for(yonNo=0;yonNo<6;yonNo++) // looking at all directions during fast marching.
% 	{
% 	  x = (int)(i - yon[yonNo][0]);
% 	  y = (int)(j - yon[yonNo][1]);
% 	  z = (int)(k - yon[yonNo][2]);
% 	  if( x>=0 && x<Size[0] && y>=0 && y<Size[1] && z>=0 && z<Size[2] && Boundary[z][y][x] == 0 && Ttag[z][y][x] == 200 ) // inside the domain and known
% 	    {
% 	      updatedDirection[0] = -yon[yonNo][0]; // from (x,y,z) to (i,j,k)
% 	      updatedDirection[1] = -yon[yonNo][1];
% 	      updatedDirection[2] = -yon[yonNo][2];
	      	      
% 	      UpdateNeighborhoodKnown( Tvalue, Ttag, F, Boundary, dx, dy, dz, Size, D, x, y, z, fm_map, trX, trY, trZ, 8, knownChangedImage, updatedDirection );
	      
	      
% 	    }
% 	}
%       time(&endKnown);
%       t_updateKnown = t_updateKnown + difftime( endKnown,startKnown );
%       //////////////////////////////////////////////
%       // ADDING NEW TRIAL POINTS //
%       time( &startTrial );
%       for(yonNo=0;yonNo<6;yonNo++)
% 	{
% 	  x = (int)(i - yon[yonNo][0]);
% 	  y = (int)(j - yon[yonNo][1]);
% 	  z = (int)(k - yon[yonNo][2]);
	  
% 	  if( x>=0 && x<Size[0] && y>=0 && y<Size[1] && z>=0 && z<Size[2] && Boundary[z][y][x] == 0 && Ttag[z][y][x] != 200 )
% 	    {
	      
% 	      UpdateNeighborhoodTrial( Tvalue, Ttag, F, Boundary, dx, dy, dz, Size, D, x, y, z, fm_map, trX, trY, trZ, 8, value_iter );
% 	    }
	  
	  
% 	}
      
%       time( &endTrial );
%       t_updateTrial = t_updateTrial + difftime( endTrial, startTrial );
%       /////////////////////////////
%     }
  
%   //CLEANING
%   for( k=0;k<Size[2];k++ )
%     {
%       for(j=0;j<Size[1];j++)
% 	{
% 	  for(i=0;i<Size[0];i++)
% 	    {
% 	      delete [] value_iter[k][j][i];
% 	      value_iter[k][j][i] = 0;
% 	    }
% 	  delete [] value_iter[k][j];
% 	  delete [] knownChangedImage[k][j];
% 	  value_iter[k][j] = 0;
% 	  knownChangedImage[k][j] = 0;
% 	}
%       delete [] value_iter[k];
%       delete [] knownChangedImage[k];
%       value_iter[k] = 0;
%       knownChangedImage[k] = 0;
%     }
%   delete [] value_iter;
%   delete [] knownChangedImage;
%   value_iter = 0;
%   knownChangedImage = 0;
  
%   for(i=0;i<8;i++)
%     {
%       delete [] trX[i];
%       delete [] trY[i];
%       delete [] trZ[i];
%       trX[i] = 0;
%       trY[i] = 0;
%       trZ[i] = 0;
%     }
%   delete [] trX;
%   delete [] trY;
%   delete [] trZ;
%   trX = 0;
%   trY = 0;
%   trZ = 0;
  
%   fm_map->chKnownX.clear();
%   fm_map->chKnownY.clear();
%   fm_map->chKnownZ.clear();
%   fm_map->trial.clear();
%   fm_map->trialC.clear();
%   delete fm_map;
%   fm_map = 0;
  
	
% }