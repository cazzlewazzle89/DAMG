
# DELETE READ SHARING DIRECTORIES OLDER THAN 30 DAYS

find /home/cwwalsh/public_html/tmp/* -type d -ctime +30 | xargs rm -rf
