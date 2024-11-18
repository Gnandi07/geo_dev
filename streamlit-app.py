import streamlit as st
import geemap.foliumap as geemap
import ee

service_account = 'xxx@ee-yyy.iam.gserviceaccount.com '
credentials = ee.ServiceAccountCredentials(service_account, './ee-yyy-json_name.json')
ee.Initialize(credentials)

st.set_page_config(layout="wide")

# Customize the sidebar
markdown = """
Web App URL: <https://earthsights.skyserve.ai>
"""

st.sidebar.title("About")
st.sidebar.info(markdown)
logo = "https://i.imgur.com/UbOXYAU.png"
st.sidebar.image(logo)

# Customize page title
st.title("SkyServe GeoAI Apps")

st.markdown(
    """
    This multi-page website hosts demo Geospatial applications by SkyServe.
    """
)

st.header("Instructions")

markdown = """
1. bla
2. blabla
"""

st.markdown(markdown)

m = geemap.Map(center=[40, -100], zoom=4)
m.add_basemap("OpenTopoMap")
m.to_streamlit(height=500)
