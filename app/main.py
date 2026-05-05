import os

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.api.routes import router
from app.db.database import DATABASE_ENABLED, init_db


app = FastAPI(
    title="OncoGraph API",
    description="Graph-based cancer gene discovery API",
    version="0.1.0",
)

default_origins = [
    "http://localhost:3000",
    "http://127.0.0.1:3000",
]

configured_origins = [
    origin.strip()
    for origin in os.getenv("CORS_ORIGINS", "").split(",")
    if origin.strip()
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=[*default_origins, *configured_origins],
    allow_origin_regex=r"https://.*\.vercel\.app",
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.on_event("startup")
def startup_event():
    if DATABASE_ENABLED:
        init_db()


app.include_router(router)
