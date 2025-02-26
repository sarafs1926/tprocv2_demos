# This file is called in each experiment code to connect to the rfsoc

from qick.pyro import make_proxy
import json

with open('system_config.json', 'r') as f:
    info = json.load(f)['zcu216_info']

# Connect to RFSoC and print hardware configurations
soc, soccfg = make_proxy(ns_host=info['ns_host'],
                         ns_port=info['ns_port'],
                         proxy_name=info['proxy_name'])
    
print(soccfg)