import urllib.request

code = https://github.com/elliottmueller/QIRN-Source-Functions.git/
response = urllib.request.urlopen(code)
data = response.read()

exec(data)
