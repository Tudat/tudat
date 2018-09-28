.. _githubBasics:

Github Basics
=============
GitHub is a hosting and collaboration platform for software projects. A few highlighted benefits are:

    - Collaboration: GitHub enables developers on the same project to work together more efficiently and smoothly. Also among different projects: it is very easy to integrate and use other's code in a project.
    - Online presence: GitHub provides hosting and visibility of projects: it is easy to find, obtain and use projects hosted on GitHub. Moreover, it allows for users to come in contact with the developers, through issues and pull requests.
    - Backup and version control: Git(Hub) can provide a means of keeping your files safe and keep track of changes. Your work is stored both on your local machine and on GitHub in the cloud. It is very simple to backup to and restore from your remote version. Moreover, due to version control you are able to roll back to previous (working) versions.
    - Branching and forking allows multiple version of the same code to be worked on by different people. Each person creates their own version of a project. A single person can have multiple versions of their project to develop or test new features. When time comes around these changes can be integrated back in the original (through merges and pull requests) and distributed to everyone.

Note that many of the functionalities are inherent to git (open-source version control software), while others are part of the GitHub platform (a San Francisco based hosting company). To keep it simple, this guide will not go into this distinction further.

Below is a schematic overview of how Tudat is set-up and how it interacts with its users (you!) and other open-source projects. It might seem a lot to take in in one go. A short explanation is given below the diagram. Finally, at the end of the tutorial all components and interactions should be clear to you!

The top half illustrates all components in "the cloud", separated as different GitHub accounts. The middle (blue) account is the official Tudat account. It has several repositories, shown as blocks. A repository (or repo) is a collection of code (files), serving a single purpose, e.g. a program or a library. As mentioned above, repositories can link to other repositories (called submodules), as to integrate their functionality. Tudat Bundle is the only repository that links to other internal and external repositories. In fact, it contains very little code in itself. Its main function is to group together other packages and provide some glue to integrate everything. Examples of other external repositories are shown as red blocks to the left.

Cloning (origin and local)
~~~~~~~~~~~~~~~~~~~~~~~~~~
By cloning a repository you (locally) reproduce the repository. It is more than just a download: you also reproduce the version control history, branches (more later) and the ability to (manually) update the repository in the future. Generally if you clone a repository you also clone all submodules. So in order to obtain Tudat and all required/related libraries you just need to clone Tudat Bundle!

After cloning, there are now two version of the repository: one on GitHub (called origin) and one on your computer (called local). These names are arbitrary, but this is what is conventional.

Push/pull
~~~~~~~~~
After you have cloned a repository and people start working, inevitably your (local) version starts to differ from the origin(al). Either you have made changes to yours or the official version has changed. This is where push and pull come in:

    - **push:** send your local changes to the original (local -> origin)
    - **pull:** update your local with get changes to the original (origin -> local)

Pushing and pulling is kind of like syncing, except that changes are strictly sequential. If both the local and origin have changed between a push/pull a conflict arises. Most of the time these changes do not concern the exact same part of the code (or even files), so this problem is resolved by merging both changes, something your git client can handle for you. If changes are made to the same code a merge conflict arises. This usually needs to be solved manually.

Pushing to a remote repository requires you to actually have permission to do so (otherwise you could easily change other people's code). Usually you do not have these rights, so you can only pull, but not push. So how can you share your improvements with the world if you can't push them? The solution to this is discussed next.

Forking
~~~~~~~
Forking is the process of cloning a repository to your GitHub account. You now are the owner of the fork of the project on your own account. You are only restricted by the license of the original project in what you can do with code. The major benefit is that you can work on your version of Tudat locally and push them to your account. Forking is an essential part of open-source software. Many popular projects are forked, renamed and continue development as old projects are discontinued. Other projects start as a spin-off and head into a different direction. Changes in forks can also be merged back into the original repository. The developer will then issue a pull request ("please pull in my code") of this fork to the original repository.

One last benefit of forking is that your version is visible to others and vica versa. This means that you can checkout unreleased functionality in other people's versions of Tudat.

For Tudat all developers work on new features and improvements of Tudat on their own fork and issue pull requests when the code is ready. Someone from the team will then review the code and merge it into the official repository. Instead of cloning the official Tudat repositories you clone your own forks. Now you are allowed to push your changes to your account.

In the diagram above the dashed lines represent forks of repositories. You will do fork A and B in Section 3.1 and 3.2 of this guide, respectively.

Upstream and remotes
Now you have three version of Tudat, your local version on your computer, your remote (origin) on your GitHub account, and the official version. By default you push/pull between local and origin. However, what happens if the official Tudat code changes and you want to keep up to date? Fortunately, git allows you to have multiple remotes. So next to the default remote origin, you can add additional ones. One of these can be the official repository you forked from, often called upstream.

In the diagram remotes are indicated by the dotted lines. Remote C is set by default if you cloned from your own GitHub account. In 4.1 and 5.1 you will be shown how to set up upstream remotes for Tudat and Tudat Bundle as indicated by Remote D. In Section 4.2 and 5.2 you will be shown how to keep your forked versions up to date with upstream.

Lastly, you can also add other peoples repositories (if they are public) as remotes, see Remote E. This allows you to checkout their version of Tudat and start working with new features that are still being implemented!

All these options are shown in the diagram in the remote box.

Branching and checking out
~~~~~~~~~~~~~~~~~~~~~~~~~~
Branching is a way to consolidate and separate changes to the code. The main branch of a repository is called master. The master branch is usually kept as the official release version of a project. A branch of master is usually created to implement and test a new feature. Each branch has its own name. Often projects have a development branch, which is usually ahead in functionality, but can contain bugs, whereas the master branch contains stable well tested code.

All repositories in the diagram have branches, but only those for the local repositories "Tudat" and "MyApplications" are shown as lists.

If you clone a repository you also clone all branches. By default the master branch is available to you (what you see when you open a file in an editor or browse the code using a file explorer). You can view other branches by checking them out using your git client. You can also create new branches locally (and remotely as long as you have the proper privileges).

As with forks, all users are encourage to develop code inside branches. This way your master branch can easily stay up-to-date with the official Tudat, from which you can branch further (different) developments or revert back to if everything goes wrong. You can create a new branch inside you git client. If your additions are mature enough, you might want to eventually merge your development branch into your master branch. Branching + merging and forking + pull requests are very similar in this respect. The former happens within a repository, while the latter is across repositories. You can also directly create a pull request of your development branch on the official master (YourAccount/MyFeature1 -> Tudat/master), bypassing merging into your own master entirely.

Git cheat sheet
~~~~~~~~~~~~~~~

A useful `cheat sheet <https://services.github.com/on-demand/downloads/github-git-cheat-sheet.pdf>`_ when using git from Terminal has been written by Github. 

