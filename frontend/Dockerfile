# Matches the version CI tests on. Node 20 cannot run the test suite —
# jsdom's undici needs an API it does not have — so a container pinned to
# it would pass here and fail the moment anyone ran the tests inside it.
FROM node:22-alpine

WORKDIR /app

COPY package.json package-lock.json* ./
RUN npm install

COPY . .

EXPOSE 5173

CMD ["npx", "vite", "--host", "0.0.0.0", "--port", "5173"]
