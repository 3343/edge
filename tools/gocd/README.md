Notification plugins require trusted SSL-certificates.
To trust the GoCD certificate in Java use the following command.
Make sure that HOST and PORT match.
HOST will also be validated as part of the auth and has to match the specs.
For example: HOST=deb, PORT=8154, if your system is called deb and GoCD's SSL runs at 8154.
```
wget https://raw.githubusercontent.com/ssbarnea/keytool-trust/master/keytool-trust
chmod +x keytool-trust
./keytool-trust HOST PORT
sudo /etc/init.d/go-server restart
```
