# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import whatapi
import secret
import time
import numpy as np

# <codecell>

def get_group_ids( search_string, api_handle=None ):
    '''
    Given a search string, returns all group ids with files containing that string.
    
    Input:
        search_string - String to search filenames for
        api_handle - what.cd api handle, default=None which means make one
    Output:
        group_ids - list of matching group ids
    '''
    if api_handle is None:
        api_handle = whatapi.WhatAPI(username=secret.username, password=secret.password)
    # Current page to list results from
    page = 1
    # List of group ids returned by search
    group_ids = []
    # Loop until out of results or error
    while(True):
        # Query what.cd API
        response = api_handle.request("browse", filelist=search_string, page=str(page) )
        # Check if there was an error
        if response['status'] != 'success':
            print "An error occurred."
            break
        # If the result set was empty (which can mean "no more results"), break
        if response['response']['results'] == []:
            break
        # Get group IDs
        for result in response['response']['results']:
            group_ids += [result['groupId']]
        # "Refrain from making more than five (5) requests every ten (10) seconds"
        time.sleep(2)
        page += 1
    return group_ids

# <codecell>

api_handle = whatapi.WhatAPI(username=secret.username, password=secret.password)

# <codecell>

acapella_ids = []
spellings = ["acapela",
             "acapella",
             "acappela",
             "accapela",
             "acappella",
             "accappela",
             "accapella",
             "accappella",
             "\"a capela\"",
             "\"a capella\"",
             "\"a cappela\"",
             "\"a cappella\""]
for spelling in spellings:
    acapella_ids += get_group_ids( spelling, api_handle )
    print spelling, len(acapella_ids)
    time.sleep(2)
np.savetxt( 'acapella_ids.txt', np.unique( acapella_ids ) )

np.savetxt( 'instrumental_ids.txt', get_group_ids( "instrumental", api_handle ) )

# <codecell>

def has_multiple_formats( group_ids, api_handle=None ):
    '''
    Given a group ID, check whether it has an instrumental in one format and an acapella in another.
    
    Input:
        group_ids - List of what.cd torrent group ids to check
        api_handle - what.cd api handle, default=None which means make one
    Output:
        valid_group_ids - list of suitable group ids
    '''
    if api_handle is None:
        api_handle = whatapi.WhatAPI(username=secret.username, password=secret.password)
    # Store group IDs which have acapella/instrumental in different formats
    valid_group_ids = []
    # Check all provided group ids
    for group_id in group_ids:
        # Query API
        time.sleep(2)
        response = api_handle.request('torrentgroup', id=group_id )
        # If there was an error, skip this group ID
        if response['status'] != 'success':
            print "An error occurred."
            continue
        # If there's only one torrent, it doesn't satisfy
        if len(response['response']['torrents']) == 1:
            continue
        # Check for vinyl format, and check whether it has instrumental or acapella (or both)
        vinyl_has_acapella = False
        vinyl_has_instrumental = False
        for torrent in response['response']['torrents']:
            if torrent['media'] == 'Vinyl':
                # Try all spellings (only need to check for part of the string)
                for spelling in ['capela', 'capella', 'cappela', 'cappella']:
                    if spelling in torrent['fileList'].lower():
                        vinyl_has_acapella = True
                # Check whether there's an instrumental version too
                if 'instrumental' in torrent['fileList'].lower():
                    vinyl_has_instrumental = True
        # If there's no vinyl format, skip it
        if (not vinyl_has_acapella) and (not vinyl_has_instrumental):
            continue
        # Check that a CD format has the opposite
        for torrent in response['response']['torrents']:
            # CD or WEB are both fine
            if torrent['media'] == 'CD' or torrent['media'] == 'WEB':
                # Try all spellings again
                for spelling in ['capela', 'capella', 'cappela', 'cappella']:
                    # Does this digital format have an acapella and the vinyl has the instrumental?
                    if spelling in torrent['fileList'].lower() and vinyl_has_instrumental:
                        valid_group_ids += [group_id]
                        print group_id
                if 'instrumental' in torrent['fileList'].lower() and vinyl_has_acapella:
                    valid_group_ids += [group_id]
                    print group_id
    # Return all unique valid group IDs
    return np.unique(valid_group_ids).tolist()

# <codecell>

both = np.unique( np.loadtxt( 'valid_group_ids.txt' ) )
for group_id in both:
    print "https://what.cd/torrents.php?id={}".format( int(group_id) )

# <codecell>

def has_both_digital( group_ids, api_handle=None ):
    '''
    Given a list of group IDs, check whether it has an instrumental AND an acapella on a digital format
    
    Input:
        group_ids - List of what.cd torrent group ids to check
        api_handle - what.cd api handle, default=None which means make one
    Output:
        valid_group_ids - list of suitable group ids
    '''
    if api_handle is None:
        api_handle = whatapi.WhatAPI(username=secret.username, password=secret.password)
    # Store group IDs which have acapella/instrumental in different formats
    valid_group_ids = []
    # Check all provided group ids
    for group_id in group_ids:
        # Query API
        time.sleep(2)
        try:
            response = api_handle.request('torrentgroup', id=group_id )
        except:
            print "An error occurred."
            continue
        # If there was an error, skip this group ID
        if response['status'] != 'success':
            print "An error occurred."
            continue
        # Check that a CD format has the opposite
        for torrent in response['response']['torrents']:
            # CD or WEB are both fine
            if torrent['media'] == 'CD' or torrent['media'] == 'WEB':
                # Want lossless torrents only
                if torrent['encoding'] == 'Lossless':
                    # If there's an instrumental
                    if 'instrumental' in torrent['fileList'].lower():
                        # Try all spellings again
                        for spelling in ['capela', 'capella', 'cappela', 'cappella']:
                            # Does this digital format have an acapella and the vinyl has the instrumental?
                            if spelling in torrent['fileList'].lower():
                                valid_group_ids += [group_id]
                                print group_id
    # Return all unique valid group IDs
    return np.unique(valid_group_ids).tolist()

# <codecell>

acapella_ids = np.array( np.loadtxt( 'acapella_ids.txt' ), dtype=np.int )
instrumental_ids = np.array( np.loadtxt( 'instrumental_ids.txt' ), dtype=np.int )
has_both = np.intersect1d( acapella_ids, instrumental_ids )

# <codecell>

has_both = has_both_digital( has_both, api_handle=api_handle )

