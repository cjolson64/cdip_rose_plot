#import requests
from bs4 import BeautifulSoup

station_url = 'http://cdip.ucsd.edu/m/deployment/station_view/'

stations_file = 'stations.html'
stations_html = open(stations_file, 'r').read()

#response = requests.get(station_url)
#stations_html = response.text

soup = BeautifulSoup(stations_html, 'html.parser')

table_html = soup.body.find('table', attrs={'id': 'StationsTable'})

print(table_html.text)
