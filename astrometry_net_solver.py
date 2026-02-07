import requests
import time
import sys

filename = r"C:\Users\jon\Desktop\ngc7008.png"
server = "http://localhost:8000"

print(filename)

# Upload
with open(filename, "rb") as f:
    r = requests.post(f"{server}/api/upload", files={"upload": f})

if r.json()["status"] == 'error':
    print("Error uploading image.")
    sys.exit(1)

job = r.json()["job_id"]

# Poll
while True:
    r = requests.get(f"{server}/api/jobs/{job}")
    status = r.json()["status"]
    if status == "solved":
        fits_url = r.json()["fits_url"]
        break
    elif status == "error":
        raise RuntimeError("Solve failed")
    time.sleep(5)

# Download solved FITS
r = requests.get(fits_url)
with open("m42_solved.fits", "wb") as f:
    f.write(r.content)
print("Solved FITS saved!")