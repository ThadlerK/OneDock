import streamlit as st
from utils import save_config

st.title("2️. Structure Preparation")

# --- POCKET SELECTION ---
st.subheader("Binding Pocket Configuration")
pocket_status = st.radio("Is the binding pocket known?", ["Known", "Unknown"], index=0)

if pocket_status == "Unknown":
    st.warning("Pocket is unknown. You need to run prediction.")
    if st.button("🔮 Go to Pocket Prediction"):
        save_config({"pocket_known": False})
        st.switch_page("pages/2_Pocket_Prediction.py")

else:
    st.info("Pocket is known. Please define coordinates.")
    c1, c2, c3 = st.columns(3)
    cx = c1.number_input("Center X", value=0.0)
    cy = c2.number_input("Center Y", value=0.0)
    cz = c3.number_input("Center Z", value=0.0)
    
    # Save coordinates instantly
    save_config({
        "pocket_known": True,
        "center_x": cx, 
        "center_y": cy, 
        "center_z": cz
    })
    
    if st.button("Confirm & Proceed to Docking ➡️"):
        st.switch_page("pages/3_Docking.py")