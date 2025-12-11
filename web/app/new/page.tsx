'use client';

import { useState } from 'react';
import { useRouter } from 'next/navigation';
import { api, ArticleCreateRequest } from '@/lib/api';
import LoadingSpinner from '@/components/LoadingSpinner';

export default function NewArticlePage() {
  const router = useRouter();
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [createdArticleId, setCreatedArticleId] = useState<string | null>(null);
  const [formData, setFormData] = useState<ArticleCreateRequest>({
    topic: '',
    target_audience: '',
    article_tone: '',
    code_ratio: '',
    theory_depth: '',
    environment: '',
  });
  const [showAdvanced, setShowAdvanced] = useState(false);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setLoading(true);
    setError(null);

    try {
      if (!formData.topic || !formData.topic.trim()) {
        setError('記事のテーマを入力してください');
        setLoading(false);
        return;
      }
      
      const article = await api.createArticle(formData);
      if (article && article.id) {
        setCreatedArticleId(article.id);
        // すぐに遷移させず、進捗ボタンで移動できるように保持
      } else {
        throw new Error('記事の作成に失敗しました: レスポンスが不正です');
      }
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : '記事の作成に失敗しました';
      setError(errorMessage);
      console.error('記事作成エラー:', err);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="max-w-3xl mx-auto">
      <h1 className="text-3xl font-bold text-gray-900 dark:text-white mb-8">新規記事作成</h1>
      
      {loading && (
        <div className="mb-4 p-4 bg-blue-100 border border-blue-400 text-blue-700 rounded">
          <div className="flex items-center">
            <LoadingSpinner />
            <span className="ml-2">記事を作成中です。しばらくお待ちください...</span>
          </div>
        </div>
      )}

      {error && (
        <div className="mb-4 p-4 bg-red-100 border border-red-400 text-red-700 rounded">
          {error}
        </div>
      )}

      {createdArticleId && (
        <div className="mb-4 p-4 bg-green-50 border border-green-200 text-green-800 rounded">
          <div className="flex items-center justify-between">
            <span>記事が作成されました。進捗を確認できます。</span>
            <button
              onClick={() => router.push(`/articles/${createdArticleId}`)}
              className="px-4 py-2 bg-green-600 text-white rounded-md hover:bg-green-700"
            >
              進捗を確認する
            </button>
          </div>
        </div>
      )}

      <form onSubmit={handleSubmit} className="bg-white dark:bg-gray-800 rounded-lg shadow p-6">
        <div className="mb-6">
          <label htmlFor="topic" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
            記事のテーマ <span className="text-red-500">*</span>
          </label>
          <textarea
            id="topic"
            value={formData.topic}
            onChange={(e) => setFormData({ ...formData, topic: e.target.value })}
            required
            rows={4}
            className="w-full px-3 py-2 border border-gray-300 rounded-md shadow-sm focus:outline-none focus:ring-blue-500 focus:border-blue-500 dark:bg-gray-700 dark:border-gray-600 dark:text-white"
            placeholder="例: PythonでQiita記事投稿を自動化するパイプライン構築"
          />
        </div>

        <div className="mb-6">
          <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
            対象読者
          </label>
          <div className="space-y-2">
            {['初心者', '中級者', '上級者', '専門家'].map((audience) => (
              <label key={audience} className="flex items-center">
                <input
                  type="radio"
                  name="target_audience"
                  value={audience}
                  checked={formData.target_audience === audience}
                  onChange={(e) => setFormData({ ...formData, target_audience: e.target.value })}
                  className="mr-2"
                />
                <span className="text-gray-700 dark:text-gray-300">{audience}</span>
              </label>
            ))}
          </div>
        </div>

        <button
          type="button"
          onClick={() => setShowAdvanced(!showAdvanced)}
          className="mb-4 text-sm text-blue-600 hover:text-blue-800 dark:text-blue-400 dark:hover:text-blue-300"
        >
          {showAdvanced ? '▼' : '▶'} オプション設定
        </button>

        {showAdvanced && (
          <div className="mb-6 space-y-4 border-t pt-4">
            <div>
              <label htmlFor="article_tone" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
                記事トーン
              </label>
              <select
                id="article_tone"
                value={formData.article_tone || ''}
                onChange={(e) => setFormData({ ...formData, article_tone: e.target.value })}
                className="w-full px-3 py-2 border border-gray-300 rounded-md shadow-sm focus:outline-none focus:ring-blue-500 focus:border-blue-500 dark:bg-gray-700 dark:border-gray-600 dark:text-white"
              >
                <option value="">選択してください</option>
                <option value="formal">フォーマル</option>
                <option value="casual">カジュアル</option>
                <option value="technical">技術的</option>
              </select>
            </div>

            <div>
              <label htmlFor="code_ratio" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
                コード比率
              </label>
              <select
                id="code_ratio"
                value={formData.code_ratio || ''}
                onChange={(e) => setFormData({ ...formData, code_ratio: e.target.value })}
                className="w-full px-3 py-2 border border-gray-300 rounded-md shadow-sm focus:outline-none focus:ring-blue-500 focus:border-blue-500 dark:bg-gray-700 dark:border-gray-600 dark:text-white"
              >
                <option value="">選択してください</option>
                <option value="low">低（10-20%）</option>
                <option value="medium">中（30-50%）</option>
                <option value="high">高（60%以上）</option>
              </select>
            </div>

            <div>
              <label htmlFor="theory_depth" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
                理論の深さ
              </label>
              <select
                id="theory_depth"
                value={formData.theory_depth || ''}
                onChange={(e) => setFormData({ ...formData, theory_depth: e.target.value })}
                className="w-full px-3 py-2 border border-gray-300 rounded-md shadow-sm focus:outline-none focus:ring-blue-500 focus:border-blue-500 dark:bg-gray-700 dark:border-gray-600 dark:text-white"
              >
                <option value="">選択してください</option>
                <option value="shallow">浅い（実践重視）</option>
                <option value="medium">中（バランス）</option>
                <option value="deep">深い（理論重視）</option>
              </select>
            </div>

            <div>
              <label htmlFor="environment" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
                想定環境
              </label>
              <input
                id="environment"
                type="text"
                value={formData.environment || ''}
                onChange={(e) => setFormData({ ...formData, environment: e.target.value })}
                className="w-full px-3 py-2 border border-gray-300 rounded-md shadow-sm focus:outline-none focus:ring-blue-500 focus:border-blue-500 dark:bg-gray-700 dark:border-gray-600 dark:text-white"
                placeholder="例: Python 3.9+, Ubuntu 22.04"
              />
            </div>
          </div>
        )}

        <div className="flex justify-end space-x-4">
          <button
            type="button"
            onClick={() => router.back()}
            className="px-4 py-2 border border-gray-300 rounded-md text-gray-700 dark:text-gray-300 hover:bg-gray-50 dark:hover:bg-gray-700"
          >
            キャンセル
          </button>
          <button
            type="submit"
            disabled={loading || !formData.topic}
            className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 disabled:opacity-50 disabled:cursor-not-allowed"
          >
            {loading ? '作成中...' : '記事を作成'}
          </button>
        </div>
      </form>
    </div>
  );
}

