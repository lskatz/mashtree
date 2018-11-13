# Contributing

Please consider contributing in the plugins directory.
If you feel like there is some new executable that you would like, then please develop it!

Ultimately please contact me through email or through issues to discuss what you are doing.
If not, then I will be pleasantly surprised with the new pull request :)

## Coding style

* Two spaces for tabs.
* Lots of comments.
* Use a main() function. I prefer only a few top-level function calls in main() rather than a lot of code.
* Have meaningful exit codes.  At the least, use 0 (zero) for a successful exit.
* All code must be compatible with the license in the Mashtree package.

## Steps for creating good issues or pull requests.

* All unit tests must pass in travis
* Create a new unit test for any new script. This is a perl script under the `t` folder.  
Unfortunately I cannot make an exception for the language here because the project is built
and tested with a perl workflow.
