from __future__ import print_function
import os.path, sys
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
from google.oauth2.credentials import Credentials
import io
from googleapiclient.http import MediaIoBaseDownload

# If modifying these scopes, delete the file token.json.
SCOPES = ['https://www.googleapis.com/auth/drive.readonly']

def main(argv):
    """Shows basic usage of the Drive v3 API.
    Prints the names and ids of the first 10 files the user has access to.
    """
    #print(argv)
    #sys.exit(1)

    imgname = argv[0]
    out_dir = argv[1]

    creds = None
    # The file token.json stores the user's access and refresh tokens, and is
    # created automatically when the authorization flow completes for the first
    # time.
    if os.path.exists('token.json'):
        creds = Credentials.from_authorized_user_file('token.json', SCOPES)
    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(
                '/home/ye6/hys_test/gdrive_test/credentials.json', SCOPES)
            creds = flow.run_console()  #run_local_server(port=43177)  #port=0
        # Save the credentials for the next run
        with open('token.json', 'w') as token:
            token.write(creds.to_json())



    service = build('drive', 'v3', credentials=creds)

    #q="mimeType = 'application/vnd.google-apps.folder'",
    # Call the Drive v3 API
    results = service.files().list(q="name = '{}'".format(imgname),
          pageSize=3,fields="nextPageToken, files(id, name)").execute()
    items = results.get('files', [])

    if not items:
        print('No files found.')
    else:
        print('Files:')
        for item in items:
            print(u'{0} ({1})'.format(item['name'], item['id']))

            if imgname==item['name']:
                request = service.files().get_media(fileId=item['id'])
                #fh = io.BytesIO()
                fh = io.FileIO(out_dir+'/'+item['name'], mode='wb')
                downloader = MediaIoBaseDownload(fh, request,chunksize=1024*1024*256)
                done = False
                while done is False:
                    status, done = downloader.next_chunk()
                    print("Downloaded %d%%." % int(status.progress() * 100),end='\r')
                
                break


if __name__ == '__main__':
    main(sys.argv[1:])