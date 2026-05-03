#!/usr/bin/env python3
    
'''
########################################################################
# AUTHOR:: AsifIqbal => AIB_EM
# USAGE :: DOWNLOAD ALL Alexandria DATABASE
########################################################################
'''

import requests
from bs4 import BeautifulSoup
from concurrent.futures import ThreadPoolExecutor
import os

base_url = "https://alexandria.icams.rub.de/data/pbe/"
save_dir = "complete_database"
os.makedirs(save_dir, exist_ok=True)

def download_file(href):
    file_url = base_url + href
    file_path = os.path.join(save_dir, href)
    try:
        with requests.get(file_url, stream=True) as r:
            r.raise_for_status()
            with open(file_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        print(f"Downloaded: {href}")
    except Exception as e:
        print(f"Failed to download {href}: {e}")

# Get the list of files
response = requests.get(base_url)
soup = BeautifulSoup(response.text, "html.parser")
file_links = [a.get('href') for a in soup.find_all('a') if a.get('href', '').endswith('.bz2')]

# Download in parallel
num_cpus = os.cpu_count()
with ThreadPoolExecutor(max_workers = num_cpus) as executor:  # Adjust max_workers as needed
    executor.map(download_file, file_links)
