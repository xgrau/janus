# HMS JANUS

## Installation

Janus is written in go. You will need to make sure that go is installed (see you system details for how to do that). 

Getting and building Janus is relatively straightforward. You can clone the repository with:

 ```bash
git clone https://git.sr.ht/~hms/janus
```

This will make a `janus` directory. You can then go into the directory and: 

```bash
mamba create -n janus conda-forge::go==1.16.5  conda-forge::nlopt==2.7.0
mamba activate janus
go build janus.go
```

This should result in downloading all the necessary libraries and compiling what is necessary, but you may need to manually do the `nlopt` compilation (see below). 

After you have ensured that `nlopt` is ok, you can go ahead and run `go build janus.go` again. Then you should be good to run `./janus`. 

Make sure that `go-nlopt` is up to date. That is the most common problem installation/run-time issue: `https://github.com/go-nlopt/nlopt`. 

Sometimes you will get an error after compiling and running an example, and this can be because of `nlopt`. 

To ensure that you have `nlopt` correctly installed, you can go to where it is installed in your system: 

```bash
cd <home>/godir/github.com/go-nlopt/nlopt
chmod 766 <nlopt folder>
bash install_nlopt.sh 
# You will need to enter your password (and if you have run this before, you may need to delete the nlopt*/build directory). 
# But after you have completed this, rerun the build of Janus.
```

For [documentation on use click here](doc/index.md).
