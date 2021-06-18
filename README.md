## Some commands
$ pip install -U jupyter-book
$ git clone https://github.com/kthpanor/echem.git
$ cd echem
$ vi docs/dft.md
$ jupyter-book build .
$ open open _build/html/index.html

## Creating final html-version
$ ghp-import -n -p -f _build/html
This requires installing ghp-import
