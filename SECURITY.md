# Security policy

## Reporting a vulnerability

Please report security problems **privately**, not as a public issue.

- Use GitHub's [private vulnerability reporting][gh] on this repository, or
- email <mohamed.merzoug.essbo@gmail.com> with `SECURITY` in the subject.

[gh]: https://github.com/Midotech31/g-synth-app/security/advisories/new

Please include what you did, what happened, and what you expected. A proof of
concept helps but is not required to report.

Expect an acknowledgement within a week. You will be told when a fix lands,
and credited in the release notes unless you would rather not be.

## What is in scope

The deployed API and the code in this repository — in particular
authentication and session handling, the per-user project boundary, request
validation, and anything that lets one account reach another's sequences.

Denial of service through a deliberately expensive request is in scope and
worth reporting: every endpoint that runs an algorithm is rate-limited and
every input field is bounded, and a way past either is a real finding.

## What is not

- Findings against a local development deployment. `config.settings.dev` has
  `DEBUG = True`, a committed secret key and permissive hosts on purpose;
  `config.settings.prod` refuses to start with any of them.
- The Ollama service behind the **Learn** page. It answers questions and
  reaches no user data, and in the supplied Docker Compose it is not exposed
  beyond the host.
- Reports that a language model gave a wrong or unhelpful answer. That is a
  quality issue, not a vulnerability — please open a normal issue.

## A note on scientific correctness

A design that is **wrong** — oligos that will not assemble, a verification
that passes when it should fail, an enzyme with the wrong cut position — is
not a security issue, but it is treated with the same seriousness. Report it
as a normal issue with the input sequence, the enzyme pair and the options
used, and it becomes a regression test.
