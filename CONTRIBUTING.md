# Code of Conduct

The following code of conduct must be followed in all places associated with this repository
itself. This includes issues, pull requests, and all files contained within the repository. The
code of conduct applies to all contributors and community members.

In addition to following GitHub's terms of service and any laws applicable to you, please:
* Do your best to properly attribute any copied, vendored, or otherwise "borrowed" code. In particular, make sure to provide LICENSE files when necessary.
* Be courteous to other contributors/community members.
* Do not make major changes to LICENSE or CONTRIBUTING.md
    * Clarifications and spelling/grammar fixes are fine; modifying the rules is not.
* Avoid vulgar or offensive language.
* Stay on-topic; all discussion, including issues and pull requests, must relate to this repository in some way, and should not focus on something other than the code and its modification.
* No malicious code/malware of any kind, including but not limited to ransomware, adware, bloatware, and spyware.

Failure to comply with this code of conduct could result in having your access to the community
restricted, such as having offending posts removed or being barred from further submissions.

You can report code of conduct violations to a maintainer (@gerudo7 or @villaa) or the
@villano-lab group via direct message or email. If you see a violation of GitHub's ToS or of a
local or federal law, please report it to the appropriate authorities first. Thanks!

# Standards

All features of `nrCascadeSim` should be tested with Travis-CI wherever possible.  This is done
through the top-level `.travis.yml` file.  You can find Travis-CI's docs
[here](https://docs.travis-ci.com/). Instructions for adding tests can be found in `nrCascadeSim/tests/README.md`.

Addition of features should be accompanied by some form of documentation outlining the use of the
new features. Similarly, if a feature is altered, documentation should be adjusted to match.

Variable and function names should relate to what the variable or function is doing.  
**Unacceptable:** `func` - Name does not describe what the function does at all.  
**Acceptable:** `calculate` - The name indicates a category (that the function calculates something), but it is still vague.  
**Best:** `integrate` - The name clearly indicates what the function does (take an integral).

*(Spot something that's not up to our standards? Submit an issue or pull request!)*

## Python Standards

Keep code neatly organized:
* Leave space before and after multi-line function definitions
* For jupyter notebooks, try to break up cells where it is sensible - for example, separating plotting cells, cells that define functions and variables, and cells that process information. 
	* In general, avoid having notebooks that consist only of a few very long cells or entirely of very many very short cells. 
* Avoid long lines of code where it's possible to break the code up into multiple lines instead.

Example:

```python
def f(a, b):
    return a*b

i = 0
if (1 == 1):
    while (n < 3):
        i += 1
        f(1,1)
```

## Documentation Standards

All functions accessible to the user should have at least one example provided, possibly more if
the usage is complex or varies significantly.

All variables and options available to the user should be clearly defined.

Documentation files should, at the end of the file, note the date corresponding to the last time
they were updated as well as the relevant version number.


# Pull Requests

Pull Requests (PRs) are created in order to submit to the owner(s) of the repository some code for
consideration. Typically that code will address some issue or improve the code in some way, we
should be clear about how we expect PRs to improve the code in our contributing documentation.
When creating the pull request you have to supply a comparison branch.  When submitting PRs,
please default to submitting to the `develop` branch and never submit directly to the `master`
branch.  This allows us to correctly manage versioning.  If you submit a pull request directly to
the `master` branch, you will be asked to change the target to `develop` (or another applicable
branch).

In order for us to be consistent with GitFlow we should adhere to the following:

**PRs submitted by the development team:** we should locally create a feature branch using the
gitflow protocol on our local machine, and then push that branch. That branch should then be
selected as the "comparison" branch for the PR. Further for the merger to be compatible with
gitflow we should define the base branch as "develop." So the steps are:

1. create a feature branch locally make some changes and push it to the remote (GitHub)
2. open a pull request with base branch `develop` and comparison branch the feature you just created `feature/XXX`
3. When you're done committing (and pushing!) to the feature branch push the button on GitHub to merge the PR back--it will merge it to develop 
4. Delete the branch on github and your local machine and add notes to the upcoming release
5. the feature will be released when the code team does the next release. 

**PR requests submitted from outside our development team:** are very similar to those from the
development team, but the team won't have access to or control over the feature branch created. It
would be created by a fork of the repository. So it looks like this:

1. create a fork of the repository with a branch dedicated to the issue (could be the local `master` we can't enforce any naming conventions there). 
2. open a PR with base branch `develop` and the comparison branch the branch on the fork you just created. 
3. When you're done committing alert the development team in the PR by using the @villaa or other tags. 
4. This will be merged back by the development team if the criteria for code improvement are met. 
5. the feature will be released when the code team does the next release. 

All PRs will be automatically by Travis-CI.  Please note whether you updated the CI or
whether no change was needed.  If for some reason a new, untested feature is implemented, but you
are unable to implement the necessary CI, explain why and how it can be manually tested.

## Release Documentation

When a PR is accepted it will be staged for a release. We will make note of all the currently
staged changes in the RELEASENOTES.md file. It is helpful, but not necessary to put a short
description under the `Next Release` section briefly describing the work done in a PR.  

## Template
The following template is not required, but if you do not use it, please be sure to include all answers to all of the questions in some other way.

**Does your pull request resolve or partially resolve an issue?** 
Yes / No.

**If Yes, which issue?** 

**Does your pull request implement code improvements?**
Yes / No.

**Does your pull request implement any breaking changes?**
Yes / No.

**If breaking changes are implemented, please describe:**

**Testing:**  
This pull request:
[ ] Alters the existing CI in some way.
[ ] Adds a new step to the CI.
[ ] Does not introduce any features that the CI could not already test.
[ ] Is not accompanied by necessary CI changes due to some limitation described below. (Please also describe how new features can be manually tested.)

See the `README.md` file in the `nrCascadeSim/tests` directory for more instructions on how to make tests. 

**Other information:**
Anything else you want to say.

# GitHub Issues

Issues fall into three categories:
* Bug report
* Feature request
* Documentation issue

When submitting issues, please be specific. We can only resolve a bug if we know what about the
program isn't working, implement a feature if we know what aspect is being improved, and clarify
documentation if we know what part is unclear.

Below are outlines for determining what your issue qualifies as and how to report them. When
submitting an issue, please specify which of these three categories you think it belongs in. We
understand that the three categories can overlap, so don't worry too much if you aren't sure if
the category you chose is appropriate. Each section also provides a template; these are just to
help people know what to write, and their use is not strictly required (although it may help us
address the issue faster).

## Bug report

When submitting a bug report, please make sure to include any information necessary for us to
reproduce the bug. If we can't reproduce it, it will be much harder to diagnose and solve the
issue.

An issue is a bug report if:
* The code does not run or only partially runs.
* The code does not build.
* A command associated with the code fails to run despite matching the documentation.
* The code takes an inordinately long amount of time to run, taking hardware into account.
* The code gives no output either to a binary file, a log file, or the terminal, when it should be giving some kind of output there.
* A command that should give consistent output gives different output each time.
* The result of a command directly contradicts what the documentation says should occur.

An issue is not a bug report if:
* The code does not interface with an environment that the documentation does not specify it will interface with. (Feature request)
* The code is missing the ability to do something you think it should be able to do, but the documentation does not specify it is capable of. (Feature request)
* The documentation is unclear but the code does not give results that directly contradict it. (Documentation issue)

### Template

I am submitting a **bug report**.

**This bug occurs in:**
make / realizeCascades / CI / specific file / etc.

**Expected behavior:**
____ should ____.

**Current behavior:**
____ instead does ____.

**Steps to reproduce:**
1. Do thing
2. Do thing
3. Result

**Other Information:**
Anything else you want to say.

**Relevant Output:**
Provide a log file, text from terminal, "No output", etc.

## Feature request

An issue is a feature request if:
* You are requesting for the code to interface in a new way that it currently does not, such as a new command or argument.
* You are proposing a particular way to increase the speed of the code.
* You are pointing out where the code could be more user-friendly.
* You are otherwise requesting for the code to do something it is not yet written to do.

An issue is not a feature request if:
* It does not affect the code, only the documentation. (Documentation issue)
* It is to fix unexpected behavior. (Bug report)
* You are providing the feature you are requesting. (Pull request)

### Template

I am submitting a **feature request**.

**The feature I am requesting is for:**
make / realizeCascades / CI / specific file / Something new / etc.

**I am requesting:**
A completely new feature / An improvement on an existing feature / etc.

**Ideas for implementation:**
(Optional)

**Other Information:**
Anything else you want to say.

## Documentation

An issue is a documentation issue if:
* A command mentioned in the documentation cannot be found.
* You would like an example and there is no similar example in the documentation.
* There is a part of the documentation you are asking for us to clarify.
* There are spelling and grammar errors, missing images, or broken links which you do not know how best to fix with a pull request.

An issue is not a documentation issue if:
* You provide fixed wording, spelling, and/or grammar for all issues you point out. (Pull request)
* The code attempts to run but fails. (Bug report)
* You are looking for a way to do something that you do not know exists and is not mentioned in the documentation. (Feature request)

### Template

I am submitting a **documentation issue**.

**The file(s) in question is/are:**
README.md / CONTRIBUTING.md / LICENSE / etc.

**The problem is in the following category/categories:**
Clarity / Examples / Broken links and images / Typos, spelling, and grammar / Undocumented Information / Out-of-date / Other

**Description of the problem:**
Describe whatever is wrong with the documentation or could otherwise be improved.

---

*Last updated 16 November, 2021, v1.0.6*
