'''Builds a custom feed from NPR podcasts 
'''

from urllib2 import urlopen
from urllib import quote
import json
import sys

def build_api_call(key, npr_id, search_string, feed_title):
    ''' builds an api call url from the provided api key, the search keywords, and the feed title
    '''
   
    url = 'http://api.npr.org/query?apiKey=' + key

    parameter1 = '&requiredAssets=audio'
    parameter2 = '&numResults=3'
    parameter3 = '&format=json'
    parameter4 = '&action=or'

    url += parameter1 + parameter2 + parameter3 + parameter4
    
    if (npr_id or search_string):
        raw_input("Hit Enter to download your podcast.")
        if npr_id:
            url += '&id=' + npr_id
        if search_string:
            url += '&searchTerm=' + quote(search_string)
        if feed_title:
            url += '&title=' + feed_title    
           
    else: 
        sys.exit("You must enter an NPR_ID, a search term, or both.")
    
    return url 
    
def call_station_api(url, feed_title, save_file):
    '''calls the NPR api with the request specified in url, returns the responds in json 
       format, and saves it under feed_title.json
    '''
    response = urlopen(url)
    j = json.load(response)
    my_feed = open(feed_title+'.json','w')
    json.dump(j, my_feed)
    my_feed.close() 
    return j
    
# NPR API key    
key = "MDE2ODA5NjY2MDE0MTE0MDgzMzdmMWNmNw001"

# Search inquiry
npr_id = raw_input("Enter comma-separated NPR IDs (visit http://www.npr.org/api/mappingCodes.php)\
 or leave blank: ")
search_string = raw_input("Enter your search string or leave blank: ")
feed_title = raw_input("What's your feed title? ")

# Build the url request
url = build_api_call(key, npr_id, search_string, feed_title)

# Executes the get request
json_obj = call_station_api(url, feed_title)

# Output the search results
print '\n\n SEARCH RESULTS \n\n'

for story in json_obj['list']['story']:
    print "TITLE: " + story['title']['$text'] + "\n"
    print "DATE: "   + story['storyDate']['$text'] + "\n"
    print "TEASER: "   + story['teaser']['$text'] + "\n"
    if 'byline' in story:
        print "BYLINE: "   + story['byline'][0]['name']['$text'] + "\n"
    if 'show' in story:
        print "PROGRAM: "   + story['show'][0]['program']['$text'] + "\n" 
    print "NPR URL: " + story['link'][0]['$text'] +'\n'           
    print "MP3 AUDIO: " +story['audio'][0]['format']['mp3'][0]['$text'] + "\n\n\n"
    
