import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
 # Title of the application
st.title("Streamlit Research Study Dashboard")
 # Create a sample data frame
data = pd.DataFrame({
 "x": np.arange(1, 101),
 "y": np.random.normal(0, 1, 100)
 })
 # Add a slider to select a subset of the data
subset_size = st.slider("Select the number of data points", 10, 100, 50)
 # Display a line plot based on the selected subset size
plt.plot(data["x"][:subset_size], data["y"][:subset_size], "-o")
st.pyplot(plt)
 # Add a text input for comments or notes
notes = st.text_area("Add your comments or notes here")