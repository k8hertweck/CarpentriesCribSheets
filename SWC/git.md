# Software Carpentry Unix Instructor Crib Sheet

## Before class:
* set up shell:
	* enlarge text size
	* `export PS1='$ '`
	* `export PROMPT_COMMAND="history 1 >> ~/Dropbox/GitHistory.txt"`
* check software installation for git AND account with GitHub
* need proxy?
* get slides setup
* remind students to register for GitHub account

## Resources
* Git for beginners: http://stackoverflow.com/questions/315911/git-for-beginners-the-definitive-practical-guide 
* fixing detached head: http://stackoverflow.com/questions/10228760/fix-a-git-detached-head 
	
## SETUP
**Objectives:** get system set up for proper attribution of your work

* Slides: run locally, Git.pdf
* open shell, can be anywhere
* specify your name, to be recorded in commits
* `git config --global user.name "k8hertweck"`
	* specify your email, to be recorded in commits
* `git config --global user.email "k8hertweck@gmail.com"`
	* specify text editor, to be used in committing
	* nano: `git config --global core.editor "nano -w"`
* `--global` means will apply for every command entered afterward
* check settings: `git config –list`
* syntax for git: `git VERB`
* find help:
	* `git config -h`
	* `git config --help`

## CREATING A REPO (WORKING LOCALLY)
**Objectives:** start tracking versions in a particular folder/directory (repository)

* start in home directory
* create directory and change directories
* `mkdir my_project`
* `cd my_project`
* tell git where to store old records and versions of file
	* `git init`
* list contents: `ls`
* list everything (including hidden files)
	* `ls -a`
* result is a subdirectory with information about project
* if this is removed, we don't have access to versioning anymore
* ask status of project
* `git status`
* output also adds helpful comments

## TRACKING CHANGES: 
**Objectives:** practice workflow (modify-add-commit), explain where information is stored

* create new file: `nano code.sh`
* add text, exit and save, `ls`, `cat mars.txt`
* show status of project: `git status`
	* draw attention to "untracked files": there's something not being tracked
* add file: `git add code.txt`
* `git status` "Changes to be committed": it's tracking, but the changes aren't recorded
* commit file: `git commit -m "creating code file"`
	* commits record to history in .git, called a revision, with short identifier
	* run without -m and will open editor to add comment (you've specified your preference earlier today)
	* good commits: brief (<50 characters), for more info, add empty line (adds in separate field)
* `git status`
* check history: `git log`
	* lists all revisions in reverse chronological order: full identifier (relate to short identifier, same initial characters),
	* when created, log message
* make another change (this is just one way): `vi code.txt`, `cat code.txt`
* check status: `git status`
* look at differences: `git diff`
	* first line: old and new version of files, similar to diff command in Unix
	* second line: labels for revisions
	* third and fourth: name of file changing
	* last lines: actual changes
* commit changes: `git commit -m "adding change"`
* but change hasn't been added!: `git add`, `git commit`
* working through staging
* new file: `vi data.txt`
* look at differences: `git diff`
* stage file: `git add data.txt`, `git diff` (no output)
* `git diff –-staged` (shows differences)
* `git commit -m`, `git status`, `git log`
* git commit workflow
* Challenge: What commands would you use to save the changes of a new file, `test.txt`, to your local Git repo?
* Challenge: Create a new Git repository on your computer called bio; Write a three-line biography for yourself in a file called me.txt, commit your changes; Modify one line, add a fourth line; Display the differences between its updated state and its original state.

## EXPLORING HISTORY
**Objectives: identify and use Git commit numbers, compare versions, restore old versions**

* prompts in git output can help you!
* compare different versions of commits
	* HEAD~1 (HEAD minus 1) and HEAD~2 refer to old commits: most recent end of chain is HEAD
	* `git diff HEAD~1 code.txt`
	* `git diff HEAD~2 code.txt`
	* `git show HEAD~2 code.txt` : includes changes and commit message
	* `git log` gives strings of digits and letters: `git diff XXXXXXXX code.txt`
	* can also just use first few characters: `git diff XXX code.txt`
* how to revert to old version?
	* make another change to code.txt, `git status`
	* `git checkout HEAD code.txt` to remove unstaged changes (default to previous committed version)
	* `git checkout XXX code.txt` to go back further in history
	* remember that you want changes before most recent commit
	* if you use `git checkout` without a file name, may end up with detached head
	* git revert
* Challenge: git checkout can be used to restore a previous commit when unstaged changes have been made, but will it also work for changes that have been staged but not committed? Make a change to code.txt, add that change, and use git checkout to see if you can remove your change.

## IGNORING THINGS
**Objectives: ignore some files, why ignoring is useful**

* sometimes backup files are created by other programs, or intermediate files we don't want to track
* create dummy files
* `mkdir results`, `touch a.dat b.dat c.dat results/a.out results/b.out`
* `git status`: lots of stuff needs to be committed, but we don't want to!
* `vi .gitignore`, add `*.dat` and `results/`
* `git status`: now .gitignore is the only thing there!
* `git add .gitignore`, `git commit -m “add ignore file”`, `git status`
* `git add a.dat`: error
* `git status --ignored`
* Challenge: Ignoring files that have already been committed retains the file in the git log. How could we remove these completely from the repo? https://help.github.com/articles/removing-sensitive-data-from-a-repository/

## REMOTES IN GITHUB
**Objectives: explain why remote repos, clone remote repos, push/pull**

* differentiating between Git and GitHub:
	* local repository: only available on your own hardware (computer)
	* remote repository: available on someone else's computer and visible to other people
* overview of GitHub
	* log in not necessary to view public repos
	* when logged in, can see timeline/newsfeed
	* starring repositories, following people
	* organizations (fredhutchio as example)
	* paid accounts can have private repos (also education discounts)
	* settings on personal page
* click icon in top right corner, create repo called `my_project`
* resulting page has info about configuring repository, it's done the `mkdir my_project`, `cd my_project`, `git init` process remotely
* connect the two repos: copy https url to clipboard
* go to local repo (in shell): `git remote add origin URL`
* `git remote -v`: check that it worked
* origin is a nickname for remote repo
* now send local changes to remote repo: `git push origin master`
* check remote repo, changes should be there (may need to refresh)
* create README file and add brief comment about the purpose of the materials, commit change, go back to terminal: `git pull origin master`
* Challenge: How is `git push` different from `git commit`?
* Challenge: You have cloned a repo owned by someone else. Can you push to and pull from that repo?

## COLLABORATING
**Objectives: collaborate pushing to common repo**

* Example repository: https://github.com/fredhutchio/guacamole
	* have default branch set to month and year for class
* issues
	* way to interact with developers about their code, request new features, report bugs
	* anyone with GitHub account can post issue
	* not always appropriate to ask questions about how to use software (see other documentation, Google groups, etc)
* Challenge: create issue
* fork repository
	* this repository is owned by fredhutchio (GitHub organization)
	* you don't have permission to work directly with the repository
	* click "fork" button in the upper right hand corner of repo
	* this creates a copy of the repository of your own that is connected with the original, but which you can edit yourself
	* relate back to branches
* edit ingredients.txt to include additional ingredients
* create new file called recipe.txt to include steps in process to making guac
* compare with original version in fredhutchio
	* can choose additional comparisons (e.g., other people, or other branches)
	* create pull request (PR): way to request that original repo will accept your changes
	* text editing options, linking people, issues, etc
	* PR etiquette: explain sufficiently the changes and why; be polite

## CONFLICTS
**Objectives: explain when conflicts occur, resolve conflicts from a merge**

* demo resolving conflict in guacamole repo
* Challenge: pair up with a partner and create pull request (PR) to send your changes to their repo
	* when comparing, your repo will be on the left and their repo will be on the right
	* don't worry yet about approving pull request

## OPEN SCIENCE
can use version control as electronic lab notebook for computational work
more open means more citation and reuse

## HOSTING
where to put code and data?
can do this yourself by purchasing domain and paying ISP to host
can also use public service
includes web interface, plus other functionalities: ability to collaborate, get DOI, academic folks can get free private repos for education 

## LICENSING
adding license and citation info
choosing licenses
licensing vs. social expectations

## WRAPPING UP
sticky notes for summary
link to materials on GitHub: http://software-carpentry.org/lessons.html 
reminders about hosting other workshops, becoming an instructor
