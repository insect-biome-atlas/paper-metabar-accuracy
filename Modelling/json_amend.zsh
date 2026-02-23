#!/bin/zsh

# Check if the input file exists
if [ ! -f "$1" ]; then
  echo "Usage: ./json_amend.zsh input_file.json"
  exit 1
fi

# Create a temporary file to store the modified JSON data
temp_file=$(mktemp)

# Initialize a variable to keep track of whether it's the first JSON object
first_json_object=true

# Add an opening square bracket at the beginning of the file
echo "[" >> "$temp_file"

# Loop through the input file, adding commas between JSON objects
while IFS= read -r line; do
  if [ "$first_json_object" = true ]; then
    # For the first JSON object, don't add a comma
    first_json_object=false
  else
    # For subsequent JSON objects, add a comma before them
    echo "," >> "$temp_file"
  fi
  # Append the JSON object to the temporary file
  echo "$line" >> "$temp_file"
done < "$1"

# Remove the trailing comma, if it exists
sed -i -e '$s/,$//' "$temp_file"

# Add a closing square bracket at the end of the file
echo "]" >> "$temp_file"

# Overwrite the input file with the modified data
mv "$temp_file" "$1"

echo "JSON file amended successfully!"

