.. _updatingTudat:

Updating Tudat
==============
During the installation as discussed in :ref:`downloadTudatBundleSmartgitDeveloper` the TudatBundle is cloned to your local system. However, this does not include the latest version of the code. This page explains how your code is updated to the newest version and how different branches of the Tudat project can be used. This tutorial assumes that you've completed the :ref:`downloadTudatBundleSmartgitDeveloper` tutorial. 

Synchronizing repositories 
~~~~~~~~~~~~~~~~~~~~~~~~~~

Step 1: Fork Tudat Bundle and Tudat repository (Do this once!)
**************************************************************
**Fork Tudat Bundle**
    Create your own copy of the Tudat Bundle under your GitHub account. Go to https://github.com/Tudat/tudatBundle and click on ``Fork``.

**Fork Tudat**
    Go to https://github.com/Tudat/tudat and click on ``Fork``.


Step 2: Synchronize with official Tudat Bundle
**********************************************
**Rename official Tudat Bundle (Do this once!)**
   In Smartgit right-click on origin choose ``Rename`` and call it "upstream". 

**Add forked Tudat Bundle as origin (Do this once!)**
   Double click on your tudatBundle repository in the ``Repositories`` tree on the left and choose ``Remote/Add...``. For the URL select https://github.com/userName/tudatBundle.git and call it "origin", click ``Add``. Right-click origin in the ``Branches`` tree on the left and select ``Pull``. A pop-up window occurs, choose ``Fetch Only``. 

**Update to the latest Tudat Bundle (Do this every time!)**
    Right-click on ``upstream`` and select ``Pull``. If a dialog pops up asking you to Rebase or Merge, select ``Rebase`` and click ``Configure``. Choose ``Fetch only``. Right-click on the ``master`` branch under the ``upstream`` entry and select ``Merge``. You will be asked how to incorporate any changes. In almost all cases ``Fast forward`` is the best option. In case there are conflicts you should do a merge. You have now succesfully updated your local Tudat Bundle to the latest version. However, you need to synchronize this change to your remote version as well.

**Synchronize with your remote (Do this every time!)**
   Right-click on ``Local Branches/master`` choose ``Push to...`` select target repository: origin (https://github.com/userName/tudatBundle.git) choose ``Push``. 

Step 3: Synchronize with official Tudat
***************************************
**Set-up official Tudat (Do this once!)**
    In the ``Repository`` window select ``tudat`` under ``tudatBundle`` and right-click the ``origin`` entry in the ``Branches`` window and select ``Properties``. Change the official remote Tudat repository url to your forked repository url: https://github/username/tudat.git. Now, all that remains is to re-add the official remote Tudat repository url. First, make sure that ``tudat`` is still selected in the ``Repositories`` window. Then, click ``Remote Add`` from the top menu. Like before fill in the official address for the remote url: ``https://github.com/Tudat/tudat.git``. Again choose "upstream" as the name.

**Update to the latest Tudat (Do this every time!)**
    Right-click on ``upstream`` and select ``Pull``. Choose ``Fetch only``. Right-click on the ``master`` branch under the ``upstream`` entry and select ``Merge``. You will be asked how to incorporate any changes. In almost all cases ``Fast forward`` is the best option. In case there are conflicts you should do a merge. You have now succesfully updated your local Tudat to the latest version. However, you need to synchronize this change to your remote version as well. It could be that you have to checkout the local ``master`` branch first if an error message pops up. Double click the local branches, ``master`` branch. Try ``Sync`` again.

Switching branches
~~~~~~~~~~~~~~~~~~
Repositories can have multiple branches. These branches are used to develop new applications or versions of the project while keeping a clear and working version of the code available in the master branch. Once the new application is verified it can be merged into the master branch for a new working version of the project. Switching between branches is simple in Smartgit:

Select a branch, for example the ``development`` branch of Tudat, select ``Check Out...``. Smartgit will ask you whether you want to create a local branch, if you intend to work on this branch and change code, create a local branch and choose ``Checkout``. If you only want to view the code inside this branch select the other option and ``Checkout``.

If you have created a local branch one can easily switch between the local branches by double-clicking on the branch and selecting ``Checkout``. 

.. warning:: Make sure to commit any changes you made in your branch before switching to a different branch. If not, you will lose all changes you've made up to the latest commit. Smartgit will likely warn you if you have made any changes to the code which have not been committed. 