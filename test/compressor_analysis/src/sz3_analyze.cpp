#define _GNU_SOURCE
#include <fenv.h>

#include "qcat_dataAnalysis.h"
#include <iostream>
#include "SZ3/api/sz.hpp"
#include "SZ3/utils/Config.hpp"
#include <stdlib.h>
#include <cmath>
#include <random>
#include <algorithm>

#define OUTFILE "/home/ac.arhammkhan/sz3_data.csv"


int main(int argc, char** argv)
{

	//turn float exception to SIGFPE, stop gdb
	feenableexcept(FE_DIVBYZERO);

	char targetFile[640];
	size_t r1=0,r2=0,r3=0,r4=0,r5=0,outSize;
	size_t nbEle = 0;
	float errorBounds[9] = {10.0, 5.0, 1.0, 1E-1, 1E-2, 1E-3, 1E-4, 1E-5, 1E-6}; //{0.01, 0.05, 0.1, 0.15, 0.2};
	int interpBlockSizes[6] = {8,16,32,64,128,256};
		
	int quantBinCnts[6];
	int n = 0;
	for(int i = 10; i <= 20; i += 2){ quantBinCnts[n] = pow(2,i); n++; }
	
	int NSAMPLES = 10;
	
	char *cmpData;
	float *decData;

	if (argc < 3) {
	
		printf("Usage: sz3_analyze [target file] [r1,r2,...]\n");
		printf("Example: sz3_analyze testfloat_8_8_128.dat 8 8 128\n");
		exit(0);
	}


	sprintf(targetFile, "%s", argv[1]);

	
	r1 = atoi(argv[2]);
	nbEle = r1;
	if(argc >= 4){
		r2 = atoi(argv[3]);
		nbEle *= r2;
	}
	if(argc >= 5){
        	r3 = atoi(argv[4]);
		nbEle *= r3;
        }

	if(argc >= 6){
		r4 = atoi(argv[5]);
		nbEle *= r4;
	}

	if(argc >= 7){
        	r5 = atoi(argv[6]);
		nbEle *= r5;
        }

	std::ifstream f(targetFile, std::ios::binary);

	printf("Dimension sizes: n5=%u, n4=%u, n3=%u, n2=%u, n1=%u\n", r5, r4, r3, r2, r1);
	//printf("%zu\n", conf.num); 

	float* data = new float[nbEle];
	f.read((char*) data, nbEle*sizeof(float));
	f.close();	

	QCAT_DataProperty* property = computeProperty(QCAT_FLOAT, data, nbEle);

	printProperty(property);	


	std::ofstream outfile(OUTFILE, std::ios::out | std::ios::in | std::ios::ate);
	
	if(!outfile.is_open())
	{
		printf("Error, can't open %s for logging", OUTFILE);
		exit(0);
	}
	
	outfile.seekp(0, std::ios::end);

	
	auto gen = std::mt19937{std::random_device{}()};

	//iterate over parameters
	//for(int i = 0; i < 20; i ++) printf("%f  ", data[i]);
	//printf("%i, %i, %i", SZ::EB_ABS, SZ::ALGO_INTERP_LORENZO, SZ::INTERP_ALGO_CUBIC);
	
	//for(int EB_MODE = SZ::EB_ABS; EB_MODE <= SZ::EB_ABS_OR_REL; EB_MODE++)
	{	
		uint8_t EB_MODE = SZ::EB_ABS;
		for(uint8_t ALG_MODE : {SZ::ALGO_LORENZO_REG, SZ::ALGO_INTERP})
		{
			for(uint8_t INTERP_MODE = SZ::INTERP_ALGO_LINEAR; INTERP_MODE <= SZ::INTERP_ALGO_CUBIC; INTERP_MODE++)
			{
				
				
				for(int k = 0; k < NSAMPLES; k++)
				{
				
					SZ::Config conf;
					if (r2 == 0) {
					    conf = SZ::Config(r1);
					} else if (r3 == 0) {
					    conf = SZ::Config(r2, r1);
					} else if (r4 == 0) {
					    conf = SZ::Config(r3, r2, r1);
					} else if (r5 == 0){
					    conf = SZ::Config(r4, r3, r2, r1);
					}
					else {
					    conf = SZ::Config(r5, r4, r3, r2, r1);
					}

					conf.errorBoundMode = EB_MODE;
					conf.cmprAlgo = ALG_MODE;
					conf.interpAlgo = INTERP_MODE;
			
			
					float eb = errorBounds[rand() % 9];
					int bs = interpBlockSizes[rand() % 6];
					int binCnt = quantBinCnts[rand() % 6];
					conf.absErrorBound = eb; //relErrorBound = eb;
					conf.interpBlockSize = bs;
					conf.quantbinCnt = binCnt;
					
					//preserve data for error calc by making a copy
                                        float* oriData = new float[nbEle];
                                        memcpy(oriData, data, nbEle*sizeof(float));
					
					outSize = 0;	
					SZ::Timer timer(true);
					cmpData = SZ_compress<float>(conf, oriData, outSize);
					double compress_time = timer.stop();

					float compress_ratio = conf.num * 1.0 * sizeof(float) / outSize;
					
					printf("outsize: %zu\t compress_ratio: %f\n", outSize, compress_ratio);


					float avg_err = 0;

					SZ::Config t_conf;
					decData = SZ_decompress<float>(t_conf, cmpData, outSize);
					int c = 0;
					float min = decData[0], max = decData[0];
					int mini = 0, maxi = 0;
					for(int i = 0; i < nbEle; i++)
					{
						
						if(!isnan(data[i]))
						{
							avg_err += abs(data[i] - decData[i]);
							if(i >= 100 && i < 110)
							{
								printf("data, decData, diff, absDiff\n%f,%f,%f,%f\n", data[i], decData[i], (data[i]-decData[i]), abs(data[i]-decData[i]));	
							}
							c++;
							if(decData[i] < min){ 
								min = decData[i];
								mini = i;
							}
							if(decData[i] > max){
								max = decData[i];
								maxi = i;
							}
						}
		
								
						
						
					}


					printf("MANUAL dec MIN/MAX: %f / %f   %i, %i\n", min, max, mini, maxi);
					
					QCAT_DataProperty* Dproperty = computeProperty(QCAT_FLOAT, decData, nbEle);
					printf("--dec data prop--\n");
					printProperty(Dproperty);

					printf("AE: %f | %zu | %i\n", avg_err, conf.num, c);
					avg_err /= conf.num;
					printf("%f\n", avg_err);


					//fname, EB_MODE, EB, CMPR_ALGO, INTERP_ALGO, interpBlockSize, quantBinCnt, size(bytes), nbEle, min, max, range, avg_value, entropy, variance, avg_err, compress ratio, compress_time  
					printf("fname, EB_MODE, EB, CMPR_ALGO, INTERP_ALGO, interpBlockSize, quantBinCnt,\ 
					size(bytes), nbEle, min, max, range, avg_value, entropy, variance, \
					Nmin, Nmax, Nrange, Navg_value, Nentropy, Nvariance, avg_err, compress_ratio, compress_time\n");

					char outStr[4096];
					sprintf(outStr, "%s, %s, %f, %s, %s, %i, %i, %i, %i, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", 
							targetFile, SZ::EB_STR[EB_MODE], eb, SZ::ALGO_STR[ALG_MODE], SZ::INTERP_ALGO_STR[INTERP_MODE], bs, binCnt,
							property->totalByteSize, nbEle, property->minValue, property->maxValue, property->valueRange, property->avgValue, property->entropy, property->zeromean_variance, 
							Dproperty->minValue, Dproperty->maxValue, Dproperty->valueRange, Dproperty->avgValue, Dproperty->entropy, Dproperty->zeromean_variance, 
							avg_err, compress_ratio, compress_time);
					
					printf("%s", outStr);
					outfile << outStr;

					delete[] oriData;
					delete[] cmpData;
					delete[] decData;
					free(Dproperty);

						
				}//samples
					
				
			}//interp

		}//cmpAlg
	}

	free(property);

	outfile.close();

	delete[] data;

	return 0;

}
