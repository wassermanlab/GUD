import os
import tempfile

import pytest

from GUD.api import app
from flask import request, jsonify

@pytest.fixture
def client():
    client = app.test_client()

    yield client


