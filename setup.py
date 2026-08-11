"""Compatibility shim for tools that still invoke ``setup.py`` directly.

Package metadata and dependencies live in ``pyproject.toml``.
"""

from setuptools import setup


setup()
