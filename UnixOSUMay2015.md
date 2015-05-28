##Software Carpentry Unix Instructor Crib Sheet
####Oklahoma State University, May 2015
####Kate L Hertweck, University of Texas at Tyler

####Before class:
* set up shell:
  * enlarge text size
  * `export PS1='$ '`
  * `export PROMPT_COMMAND="history 1 >> ~/Dropbox/UnixHistory.txt"`
* students check software installation: Unix, git (with XCode) and Python for Day 1

**Intro to Software Carpentry:** http://swcarpentry.github.io/slideshows/introducing-software-carpentry/index.html#slide-0 

####Checklist for class:
* class website: http://ouinformatics.github.io/2015-05-27-osu 
* etherpad: https://etherpad.mozilla.org/2015-05-27-osu
* Kate's history: https://www.dropbox.com/s/femrf7if1swjno6/UnixHistory.txt?dl=0
* review of pre-class survey
  * mixed audience, both novice/experienced and discplines
  * if material is review for you, help by keeping notes on etherpad and helping your neighbor 
  * scientific computing has a long history of being self taught, so most instructors even learn something new

####SETUP

* slides: http://swcarpentry.github.io/shell-novice/motivation.html
* most tasks in the shell can be done with mouse on Desktop. Why do anything differently?
* motivation
  * Unix = old school system
  * combining powerful tools together using minimal keystrokes
  * automating repetitive tasks: moving files
  * required to use high performance computing systems
* terms: file, directory/folder
* get data
  * http://swcarpentry.github.io/shell-novice/shell-novice-data.zip, move to Desktop, double click to unzip (if not already done), folder named “data”, can get there using: cd && cd Desktop/shell-novice/data
  * cd, git clone https://github.com/ouinformatics/osu-data 

####INTRODUCING THE SHELL
**Objectives:** orient to shell and how it relates to the computer, understand the benefit of CLI

* what computers do:
  * run programs
  * store data
  * communicate with each other
  * communicate with us → today you'll learn a new way of doing this
* terms:
  * command line interface: CLI
  * graphical user interface: GUI
* how it works:
  * you type something
  * computer reads it
  * executes command
  * prints output
  * use command shell to make this happen: this is the interface between user and computer
* bash: Bourne again shell, most commonly used, default on most modern implementations
* example data:
  * our friend Nelle has six months worth of survey data collected from the North Pacific
  * 300 samples of goo
  * pipeline:
    * assay abundance of 300 proteins, each sample has one output file with one line for each protein
    * calculate statistics for each protein separately using program goostat
    * compare statistics for proteins using program called goodiff
    * write up results and submit by end of month
  * if enters all commands by hand, will need to do 45,150 times. 
  * What can she do instead?

####FILES AND DIRECTORIES
**Objectives:** paths, learn basic commands for working with files and directories, learn syntax of commands, tab-completion

* prompt: `$` indicates computer is ready to accept commands
* `whoami` press enter
  * finds program
  * runs program
  * displays program's output
  * displays new prompt
* `pwd` print working directory, in this case it is also the home directory
* root directory: holds everything else, begins with slash `/`
* structure of directories: nested
* slashes can also be a separator between names
* `ls` listing, prints names of files and directories in current directory and prints in alphabetical order
* file suffixes are informative and necessary for computers (and us!) to interpret them
* `ls -F` adds trailing / to names of directories
  * spaces and capitalization in commands are important!
  * `-F` is an option, argument, or flag
* `ls -F data` to list files within data
* `ls -F /data` what is the difference here? this command will search the same directory in root, regardless of what working directory is
* `ls` again folders in current directory
* `cd data` change directory to data folder
* `pwd`
* `ls -F` 
* `cd ..`  go up one level in file hierarchy 
  * `..` is special directory
  * can also use absolute paths
* `ls -F -a` to see hidden files, including `.` and `..`
* file organization: 
  * `ls north-pacific-gyre/2012-07-03/`
  * tab completion
* Challenge questions 1 and 2 through Socrative
  * teacher room: PAW5AYWM
  * start quiz, teacher paced, disable student names

####CREATING THINGS
**Objectives:** create directory hierarchy that matches given diagram, create files, look in folders, delete folders

* go back to Nelle's home directory (how?)
* `pwd`, `ls -F`
* `mkdir thesis`
* `ls -F`
* `ls -F thesis`
* move into thesis
* `nano draft.txt`
  * creates file, opens text editor
  * write text
  * use `Control+O` to save file shorthand is `^O`)
  * `Control+X` to exit
* `ls`, `rm draft.txt`, `ls`
  * where does file go?
* recreate file then move up one directory
* `nano draft.txt`
* `cd ..`
* removing a directory:
  * `rm thesis`: error
  * `rmdir thesis`: still get error
  * `rm thesis/draft.txt`
  * `rmdir thesis`
  * could've also used `rm -r thesis`, but that can be dangerous!
* `mkdir thesis`: recreate
* `nano thesis/draft.txt`
* `ls thesis`, `mv thesis/draft.txt thesis/quotes.txt`, `ls thesis`
  * note: `mv` works on directories as well
* `mv thesis/quotes.txt .`
* `ls thesis`
* `ls quotes.txt`
* `cp quotes.txt thesis/quotations.txt`
* `ls quotes.txt thesis/quotations.txt`
* `rm quotes.txt`, `ls quotes.txt thesis/quotations.txt`
* Socrative questions 3 and 4 

####PIPES AND FILTERS
**Objectives:** redirect command output to file, construct pipelines

* `ls molecules`, `cd molecules` 
* `wc *.pdb`
  * word count: lines, words, characters
  * `*` is a wildcard, it matches anything (zero or more characters, there are others)
* `wc -l *.pdb` only report number of lines
* `wc -l *.pdb > lengths.txt : send output to new file named lengths.txt
* `>` will overwrite previous file
  * what does `>>` do?
* can't remember how wc reports? use `man wc` (`q` to exit), `wc -h`, or `wc –help` (this should work for most unix commands), also google man wc
* ls lengths.txt, cat lengths.txt (for concatenate)
* `sort -n lengths.txt` sort by first column, using numerical order
* does not change file, just prints output to screen
* `sort -n lengths.txt > sorted-lengths.txt`
* arrow up to recall last command
* `head -1 sorted-lengths.txt` final result: which one file is shortest?
* `sort -n lengths.txt | head -l` : vertical bar is a pipe, which sends output of command on left as input to command on right
  * `head` prints specified number of lines from top of file
* `wc -l *.pdb | sort -n | head -1` programming model: pipes and filters
  * note: you only enter the original files once!
  * standard input, standard output
* nelle's pipeline:
  * start in her home directory (`users/Nelle`)
  * `cd north-pacific-gyre/2012-07-03`
  * all files should contain same amount of data
  * any files contain too little data?
  * `wc -l *.txt | sort -n | head -5`
  * any files contain too much data?
* `wc -l *.txt | sort -n | tail -5`
* file marked with Z? outside naming convention, may contain missing data
  * `ls *Z.txt`
* records note no depth recorded for these samples
* may not want to remove, but will later select all other files using *[AB].txt
* Socrative questions 5 and 6

####LOOPS
**Objectives:** write loops that apply commands to series of files, trace values in loops, explain variables vs values, why spaces and punctuation shouldn't be used in file names, history, executing commands again

* what if you wanted to perform the same commands over and over again on multiple files?
* go to creatures directory `users/nelle/creatures`
* may try: `mv *.dat original-*.dat` bur doesn't work
* you can perform these operations using a loop
* example looking at first three lines in each file
```
for filename in basilisk.dat unicorn.dat
  do 
    head -3 $filename
  done
```
  * what does this look like when you arrow back up?
  * explain syntax: filename is variable, what does it stand for? how is it represented later?
  * shell prompt changes, if you get stuck, use control+C to get out
  * can specify whatever variable name you want
  * why might it be problematic to have filenames with spaces?
* you can include multiple commands in a loop:
```
for filename in *.dat
do
echo $filename
head -100 $filename | tail -20 
done
```
  * use of wildcard. what does echo do? why is this useful for loops?
* write a for loop to resolve the original problem of creating a backup (copy of original data)
  *strategy: `echo` command before running final, to make sure loop is functioning the way you expect
* nelle's example:
  * `cd north-pacific-gyre/2012-07-03`
  * check: `for datafile in *[AB].txt; do; echo $datafile; done`
  * add command: `for datafile in *[AB].txt; do; echo $datafile stats-$datafile; done`
  * add command: `for datafile in *[AB].txt; do; goostats $datafile stats-$datafile; done` (kill job using ^C)
  * add echo: `for datafile in *[AB].txt; do; echo $datafile; goostats $datafile stats-$datafile; done`
  * tab completion: move to start of line using `^A` and end of line `^E` (option with arrows to move by one word)
* `history`: see old commands, find line number (repeat using !number)
* Socrative questions 7 and 8

####SHELL SCRIPTS
**Objectives:** write shell script to run command or series of commands for fixed set of files, run shell script from command line, write shell script to operate on set of files defined on command line, create pipelines including user-written shell scripts 

* go back to molecules in nelle's directory
* create file called `middle.sh` and add this command: `head -15 octane.pdb | tail -5`
* `bash middle.sh`
  * `.sh` means it's a shell script 
  * very important to make these in a text editor, rather than in word!
* edit `middle.sh` and replace file name with `"$1"` 
  * quotations accommodates spaces in filenames
* `bash middle.sh octane.pdb`, should get same output
  * try another file: `bash middle.sh pentane.pdb`
* edit `middle.sh` with `head “$2” “$1” | tail “$3”`
* `bash middle.sh pentane.pdb -20 -5`
* to remember what you've done, and allow for other people to use: add comments to top of file
```
#select lines from middle of a file
#usage: middle.sh filename -end_line -num_lines
```
  * explain comments
* what if we wanted to operate on many files? create new file: `sorted.sh`
* `wc -l “$@” | sort -n`
* `bash sorted.sh *.pdb ../creatures/*.dat`
* add comment!
* save last few lines of history to file to remember how to do work again later: 
  * `history | tail -4 > redo-figure.sh`
  * `history | tail -5 | colrm 1 7` (1-7 characters)
* nelle problem
  * run goostats on all data files
  * `do-stats.sh`:
```
#calculate reduced stats for data files at J = 100 C/bp
  for datafile in “$@”
    do 
      echo $datafile
      bash goostats -J 100 -r $datafile stats-$datafile
  done
```
  *   `bash do-stats.sh *[AB].txt`
  * or just report `bash do-stats.sh *[AB].txt | wc -l`
* Socrative question 9

####FINDING THINGS
**Objectives:** grep to select lines in text which match patterns, find to find files whose names match patterns, nesting files, text vs binary files

* move to writing subdirectory
* `cat haiku.txt`
* `grep not haiku.txt` : find lines that contain “not”
* `grep day haiku.txt` : find lines that contain “day”
* `grep -w day haiku.txt` : searches only for whole words
* `grep -n it haiku.txt` : includes the numbers on lines that match
* `grep -n -w the haiku.txt` : combine flags
* `grep -n -w -i the haiku.txt` : make case insensitive
* `grep -n -w -v the haiku.txt` : invert selection, only lines that do NOT contain the
* the real strength of grep, and the origin of its name, is “regular expressions,” which describes ways of programmatically describing text strings/search patterns, but we don't have time to cover these today (many awesome lessons and tutorials online)
* difference between grep and find
* `find . -type d` : look for things that are directories in given path
* `find . -type f` : look for files instead
* find is automatically recursive (keeps drilling down into file hierarchy)
  * can specify depth: `find . -maxdepth 1 -type f`  
  * or `-mindepth`
* can match by name: `find . -name *.txt` (will only give one filename! expands name prior to running command)
* correct way: `find . -name '*.txt'`
* find similar to list, but has more refined parameter searching
* can combine together: count all lines in a group of files
  * `wc -l $(find . -name '*.txt')`
  * nesting or subshell
* equivalent command: `wc -l ./data/one.txt ./data/two.txt ./haiku.txt`
* can also combine find and grep: find .pbd files that contain Iron
  * `grep FE $(find .. -name '*.pdb')`
* today we've only talked about text files, what about images, databases, etc? those are binary (machine readable)
* Socrative question 10

####END CLASS
* stop shell script output
