import nox
import os

nox.options.sessions = ["black", "isort", "flake8", "tests-3.11"]
nox.options.needs_version = ">= 2024.3.2"


@nox.session
def flake8(session):
    session.install("flake8")
    session.run(
        "flake8", "--max-line-length=160", "--exclude", ".nox,*.egg,build,data",
        "--select", "E,W,F", "."
    )


@nox.session(tags=["style", "fix"])
def black(session):
    session.install("black")
    session.run("black", "vexo")


@nox.session(tags=["style", "fix"])
def isort(session):
    session.install("isort")
    session.run("isort", "vexo")


@nox.session()
def mypy(session):
    session.env['PYTHONPATH'] = os.getcwd()
    session.install("mypy")
    session.install("-r", "requirements.txt")
    session.run("mypy", "--ignore-missing-imports", "--no-namespace-packages", "--exclude=*site-packages*", "vexo")


@nox.session(python=["3.8", "3.9", "3.10", "3.11", "3.12"])
def tests(session):
    session.env['PYTHONPATH'] = os.getcwd()
    session.install("-r", "requirements.txt")
    session.run("py.test", "vexo/", *session.posargs)
