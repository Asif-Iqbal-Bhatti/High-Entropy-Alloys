#!/usr/bin/env python
	
from tqdm import tqdm
import multiprocessing as mp
from joblib import Parallel, delayed
import pickle, os, subprocess, re, sys
import multiprocessing as mp
from multiprocessing import Process

exe = 'mcsqs'	
num_cores = 40
print(num_cores)

def run_mcsqs(k):
	subprocess.Popen(exe + ' -rc -ip ='+ str(k), stdout=subprocess.PIPE, shell=True).communicate()
	#subprocess.run(exe + ' -rc -ip ='+ str(k), check=True, capture_output=True, shell=True)

if __name__ == "__main__":


	Parallel(n_jobs=-1, prefer="threads", batch_size='auto')(delayed(run_mcsqs)(k) for k in tqdm(range(num_cores)))
	

