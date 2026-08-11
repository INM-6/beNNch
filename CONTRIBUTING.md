# Contributing to beNNch

Thanks for your interest in beNNch. This document explains how to report
problems, propose changes, and what to expect in return.

## Project status

beNNch is actively developed. Issues and pull requests are welcome and will be reviewed. Bug reports and corrections to the documentation are particularly welcome.

If you are building on beNNch for your own infrastructure, you are very welcome to fork it. Please keep the license and attribution intact (see [Licensing](#licensing)).

## Which version are you working on?

The `main` branch contains **version 2.0 (CI-beNNch)** — the continuous benchmarking system built around GitLab CI.

The original beNNch benchmarking workflow is preserved at the tags [`v1.0`](https://github.com/INM-6/beNNch/tree/v1.0) and [`v1.1`](https://github.com/INM-6/beNNch/tree/v1.1). If your issue concerns 1.x behaviour, say so explicitly and name the tag, otherwise we will assume you mean 2.0.

## Reporting a bug

Open an issue at https://github.com/INM-6/beNNch/issues and include:

- the version or tag (`v2.0`, `v1.1`, or a commit SHA);
- the machine configuration you used (one of `.gitlab-ci/machine-configs/`, or a description of your own);
- the benchmark configuration YAML, with any credentials or private paths removed;
- what you expected to happen and what actually happened;
- the relevant part of the CI job log or the JUBE output.

Because the benchmark jobs run against a private GitLab instance and require a JSC compute-time allocation, we usually cannot reproduce a failure on your infrastructure. The more of the above you provide, the more likely we can help.

## Proposing a change

1. Fork the repository and create a branch off `main`.
2. Keep the change focused — one logical change per pull request.
3. Follow the conventions already in the file you are editing.
4. Every new source file must carry the GPL header and SPDX identifier used
   throughout the repository:

   ```
   SPDX-License-Identifier: GPL-3.0-or-later
   ```

5. Python code targets Python 3 and follows PEP 8.
6. Open a pull request against `main` describing what you changed and why, and how you verified it.

### Verifying a change

There is no automated test suite in this repository; the pipeline definitions are exercised by running them. If you change anything under
`.gitlab-ci/scripts/`, state in the pull request how you checked it. For example, the generated pipeline YAML before and after, or a benchmark run that completed on a specific machine.

## Citing and credit

If you use beNNch in published work, please cite it (see [`CITATION.cff`](CITATION.cff)).

Substantial contributions are acknowledged in `CITATION.cff` and `codemeta.json`. If you contribute something substantial and would like to be listed, say so in your pull request and include your ORCID.

## Licensing

beNNch is licensed under the GNU General Public License v3.0 or later (`GPL-3.0-or-later`) — see [LICENSE](LICENSE). By contributing, you agree that your contribution is licensed under the same terms.

## Contact

For questions that do not fit an issue, contact the maintainer listed in [`CITATION.cff`](CITATION.cff).
