FROM python:3.8-slim-buster

# Set the working directory in the container to /app
WORKDIR /app

# Add the current directory contents into the container at /app
ADD . /app

#RUN pip install --no-cache-dir -r requirements.txt

# Run factorize.py when the container launches
CMD ["python", "calculus.py"]