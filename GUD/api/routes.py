from GUD.api import app

@app.route('/')
def index():
    return 'HOME'